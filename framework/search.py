from .sdp import SyndromeDecoding
from .sdith_hypercube import HypercubeSDitH
from .sdith_threshold import ThresholdSDitH
from math import floor, ceil, log2
from math import comb as binomial

class Search:

    @staticmethod
    def run(with_sss=False, **kwargs):
        """ Perform an exhaustive search to find
            the shortest signature size with the given
            constraints.

            Return: a couple (size, variant)
                where "size" is the best achieved size (in bytes)
                and "variant" is a SDitH object which describes
                    the best parameter sets.

            Mandatory parameters:
              - the code length: "n"
              - the code dimension: "k"
              - the mimimal ISD cost: "lda"
              - the forgery cost: "kappa"

            Optional parameters:
              - the field size |F_sd|: "q"
                    by default,  q=256
              - the code weight: "w"
                    by default: w=-1
                    [Remark] if the value is negative, then it means that it is relative
                      to the GV bound. For example, "-2" means "GV-2".
              - the extension degree between F_sd and F_poly: "ext1"
                    by default: the minimal value such that |F_poly| >= m
              - the extension degree between F_poly and F_point: "ext2"
                    by default: the minimal value such that |F_points| >= 2^24
              - the number of evaluations: "t"
                    by default: t=1
              - the number of parties: "N"
                    by default: N=256
              - the number of iteration: "tau"
                    by default: the minimal value such that the forgery cost is above 2^kappa.
              - a score function: "get_score". By default, it is the zero function.
                    if two parameter sets lead to the same signature size, it applies the score
                    function, and keep the parameter set with the higher score.

            Order of the selection: q, n, k, w, ext1, ext2, t, N, tau

            Format of the parameters:
              - it can be an integer: the value is fixed
              - it can be a list: the search will try all the values in the list
              - it can be a function which returns a list:
                  the function will receive a dictionary with the already-selected parameters 
                  and must return a list. The search will then try all the values in this list.
                  It enables to have a dynamic list (a list which depends on the previous parameters).

            Example:

                size, variant = Search.run(
                    kappa=128,
                    lda=148,
                    q=256,
                    n=range(230,250),
                    k=range(115,140),
                    w=range(-10,0),
                    N=256,
                    t=5
                )

                will search the best parameter choice for
                    q=256, n between 230 and 250, k between 115 and 140
                    w between GV-10 and GV, N=256, t=5, kappa=128, lda=148
                    (and "ext1" and "tau" optimal)
        """

        estimate_peters_isd = kwargs.pop('estimate_peters_isd', True)
        estimate_lee_brickell_isd = kwargs.pop('estimate_lee_brickell_isd', True)

        def aux(lst, params):
            if len(lst) == 0:
                # When all the parameters are selected,
                #   estimate the signature size
                N = params['N']
                tau = params['tau']
                variant = params['variant']
                if with_sss:
                    variant.set_tradeoff(N,tau, params['ell'])
                else:
                    variant.set_tradeoff(N,tau)
                size = variant.get_sig_size()[1] # Take the average
                return size, variant

            # When it remains at least one parameter to select
            key, options = lst[0]
            lst = lst[1:]

            # Build the list of possible values for the current parameter
            is_list = True
            try:
                # Test if it is a number
                int(options)
                options = [options]
            except:
                try:
                    # Test if it is a (static) list
                    options = list(options)
                except:
                    try:
                        # Test if it is a dynamic list (via function)
                        options = list(options(params))
                    except:
                        is_list = False
            
            if is_list:
                best_size = None
                best_variant = None
                for value in options:
                    new_params = params.copy()
                    if key == 'w':
                        # If 'w' is negative, scale according to GV
                        if value <= 0:
                            value += new_params['gv']
                    new_params[key] = value
                    if key == 'k':
                        # After choosing q, n and k, let directly compute the GV distance
                        new_params['gv'] = floor(SyndromeDecoding.compute_max_weigth_for_target(
                            params['q'],
                            params['n'],
                            new_params['k'],
                            params['nb_additional']
                        ))
                    elif key == 'w':
                        # After choosing the SD instance, let compute the ISD cost
                        #   and abort when it is too small.
                        try:
                            new_params['sd'] = SyndromeDecoding(
                                params['q'],
                                params['n'],
                                params['k'],
                                new_params['w'],
                                d=params['d']
                            )
                        except AssertionError:
                            continue
                        cost1 = new_params['sd'].get_cost_peters_isd() if estimate_peters_isd else params['lda']
                        cost2 = new_params['sd'].get_cost_lee_brickell_isd() if estimate_lee_brickell_isd else params['lda']
                        if min(cost1, cost2) < params['lda']:
                            continue
                    elif key == 't':
                        # After choosing the parameter about MPC protocol, let compute the
                        #   the false positive rate
                        if with_sss:
                            new_params['variant'] = ThresholdSDitH(
                                params['sd'],
                                new_params['t'],
                                params['ext1'],
                                params['ext2'],
                                kappa=params['kappa']
                            )
                        else:
                            new_params['variant'] = HypercubeSDitH(
                                params['sd'],
                                new_params['t'],
                                params['ext1'],
                                params['ext2'],
                                kappa=params['kappa']
                            )
                        new_params['variant'].get_false_positive_probability() # load in cache
                    size, variant = aux(lst, new_params)
                    if size is None:
                        continue
                    if (best_size is None) or (best_size > size):
                        best_size = size
                        best_variant = variant
                    elif (best_size == size):
                        previous_score = get_score(best_variant)
                        current_score = get_score(variant)
                        if current_score > previous_score:
                            best_size = size
                            best_variant = variant
                return best_size, best_variant
            
            else:
                # Default strategy for the current parameter
                assert options is None, (key, options)

                if key == 'ext1':
                    q = params['q']
                    n = params['n']/params['d']
                    ext1 = 1
                    while q**ext1 < n:
                        ext1 += 1
                    new_params = params.copy()
                    new_params[key] = ext1
                    return aux(lst, new_params)

                elif key == 'ext2':
                    fpoly = params['q']**params['ext1']
                    ext2 = 1
                    while fpoly**ext2 < 2**24:
                        ext2 += 1
                    new_params = params.copy()
                    new_params[key] = ext2
                    return aux(lst, new_params)

                elif key == 'tau':
                    N = params['N']
                    if with_sss:
                        ell = params['ell']
                        tau = ceil(params['kappa'] / log2(binomial(N,ell)))
                        params['variant'].set_tradeoff(N,tau,ell)
                    else:
                        tau = ceil(params['kappa'] / log2(N))
                        params['variant'].set_tradeoff(N,tau)
                    while params['variant'].get_signature_security() < params['kappa']:
                        tau += 1
                        if with_sss:
                            params['variant'].set_tradeoff(N,tau,params['ell'])
                        else:
                            params['variant'].set_tradeoff(N,tau)
                    new_params = params.copy()
                    new_params[key] = tau
                    return aux(lst, new_params)

                else:
                    raise NotImplementedError('No default rule for {}'.format(key))
                
        # The list of parameter selection
        kappa = kwargs.pop('kappa')
        lda = kwargs.pop('lda')
        d = kwargs.pop('d', 1)
        nb_additional = kwargs.pop('nb_additional', 1)
        get_score = kwargs.pop('get_score', lambda x: 0)
        if with_sss:
            lst = [
                ('q', kwargs.pop('q', 256)),
                ('n', kwargs.pop('n')), # No default
                ('k', kwargs.pop('k')), # No default
                ('w', kwargs.pop('w', [-1])),
                ('ext1', kwargs.pop('ext1', None)),
                ('ext2', kwargs.pop('ext2', None)),
                ('t', kwargs.pop('t', 1)),
                ('N', kwargs.pop('N', 256)),
                ('ell', kwargs.pop('ell', 1)),
                ('tau', kwargs.pop('tau', None)),
            ]
        else:
            lst = [
                ('q', kwargs.pop('q', 256)),
                ('n', kwargs.pop('n')), # No default
                ('k', kwargs.pop('k')), # No default
                ('w', kwargs.pop('w', [-1])),
                ('ext1', kwargs.pop('ext1', None)),
                ('ext2', kwargs.pop('ext2', None)),
                ('t', kwargs.pop('t', 1)),
                ('N', kwargs.pop('N', 256)),
                ('tau', kwargs.pop('tau', None)),
            ]
        assert len(kwargs) == 0, 'Unknown parameters: {}'.format(list(kwargs.keys()))

        # Launch the exhaustive search
        return aux(lst, {'kappa': kappa, 'lda': lda, 'd': d, 'nb_additional': nb_additional, 'get_score': get_score})
