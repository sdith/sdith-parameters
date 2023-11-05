from .sdith_hypercube import HypercubeSDitH
from math import log2, ceil, floor
from math import comb as binomial

class ThresholdSDitH(HypercubeSDitH):
    def __init__(self, *args, **kwargs):
        self.ell = None
        super().__init__(*args, **kwargs)

    def set_tradeoff(self, N, tau, ell):
        assert N <= self.sd.q
        self.N = N
        self.tau = tau
        self.ell = ell

    def get_parameters(self, as_tuple=False):
        res = super().get_parameters(as_tuple=as_tuple)
        if as_tuple:
            return tuple(list(res) + [self.ell])
        else:
            res['ell'] = self.ell
            return res


    def get_false_positive_probability(self):
        if self.p is None:
            # Parameters
            (q, n, _, w, d, t, ext1, ext2, _, _, _) = self.get_parameters(as_tuple=True)
            self.p = self._compute_false_positive_probability(q, n, w, d, t, ext1, ext2)
        return self.p

    def can_use_same_unif(self):
        (_, _, _, _, _, _, _, _, N, tau, ell) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        p_ = tau*p*binomial(N, ell+1) # Conservative
        return (log2(p_)) <= -self.kappa

    def get_sig_size(self):
        """ Return the signature size in bytes. More precisely, it returns (maxi, avg, std)
            where "maxi" is the maximum size, "avg" is the mean size, "std" is the standard
            deviation.
        """
        (q, _, k, w, d, t, ext1, ext2, N, tau, ell) = self.get_parameters(as_tuple=True)
        same_unif = self.can_use_same_unif()

        # Components
        dig = 2*self.kappa # Digest
        salt = 2*self.kappa # Salt
        lq = ceil(log2(q))

        # Subparts
        plaintext_size = k*lq
        poly_size = w*lq*ext1
        bn20_uni_cost = t*lq*ext1*ext2

        inputs = plaintext_size + 2*poly_size + bn20_uni_cost
        comm = 2*d*bn20_uni_cost
        unif = 2*d*bn20_uni_cost

        # Upper bound on the signature size
        nb_max_open_leaves = min(
            BinaryTree.get_upper_bound_num1(N-ell,N),
            BinaryTree.get_upper_bound_num2(N-ell,N)
        )
        bitsize_maxi = dig + salt + tau*(
            dig*nb_max_open_leaves + ell*(inputs+unif)
        )
        if same_unif:
            bitsize_maxi += comm
        else:
            bitsize_maxi += tau*comm
        bitsize_maxi += 32 # Signature size (encoded on uint32_t)
        size_maxi = ceil(bitsize_maxi/8)

        # Mean of the signature size
        bitsize_avg = dig + salt + tau*(
            dig*BinaryTree.get_nb_leaves(N-ell,N)
            + ell*(inputs+unif)
        )
        if same_unif:
            bitsize_avg += comm
        else:
            bitsize_avg += tau*comm
        bitsize_avg += 32 # Signature size (encoded on uint32_t)
        size_avg = ceil(bitsize_avg/8)

        # Standard deviation of the signature size
        #   -> Not Implemented
        size_std = -1

        return size_maxi, size_avg, size_std

    def get_signature_security(self):
        """ Return the security of the signature in bits """
        (_, _, _, _, _, _, _, _, N, tau, ell) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        p *= binomial(N, ell+1) # Conservative
        return self._compute_forgery_cost(p, binomial(N, ell), tau)

    def print(self, in_bytes=False, new_line=True, with_sd_hardness=False):
        text = []
        (q, m, k, w, d, t, ext1, ext2, N, tau, ell) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        nb_solutions = self.sd.get_nb_solutions()
        same_unif = self.can_use_same_unif()

        # SD instance
        if d == 1:
            text.append('SD=({}, {}, {}, {}), nb=1+{:.4f}'.format(q,m,k,w,nb_solutions-1))
        else:
            text.append('SD=({}, {}, {}, {}), {}-split, nb=1+{:.4f}'.format(q,m,k,w,d,nb_solutions-1))

        # MPC Protocol parameters
        text.append('MPC=(t={},ext1={},ext2={})'.format(t,ext1,ext2))

        # Trade-off
        text.append('Tradeoff=(N={},tau={}, ell={}), Mono={}, kappa={}'.format(N,tau,ell,same_unif,self.kappa))

        # False positive probability
        text.append(' - False positive probability: 2^{:.1f}'.format(log2(p)))

        # Security
        text.append(' - Signature security: {:.2f} bits'.format(
            self.get_signature_security()
        ))
        if with_sd_hardness:
            cost_peters = self.sd.get_cost_peters_isd()
            cost_lee_brickell = self.sd.get_cost_lee_brickell_isd()
            text.append(' - SD security: {:.2f} bits (Peters: {:.2f}, Lee-Brickell: {:.2f})'.format(
                min(cost_peters,cost_lee_brickell),
                cost_peters,
                cost_lee_brickell,
            ))

        # Sizes
        size_maxi, size_avg, size_std = self.get_sig_size()
        if in_bytes:
            text.append(' - Size: {} B (std={},max={})'.format(
                round(size_avg), round(size_std), round(size_maxi)
            ))
        else:
            factor = 1024
            text.append(' - Size: {:.2f} KB (std={:.2f},max={:.2f})'.format(
                size_avg/factor, size_std/factor, size_maxi/factor
            ))

        # Print
        print('\n'.join(text))
        if new_line:
            print()


class BinaryTree:
    @staticmethod
    def get_upper_bound_num1(nb_revealed, nb_committed):
        """ Get Upper Bound Num 1:
            Good bound for 'nb_revealed' small """
        return nb_revealed

    @staticmethod
    def get_upper_bound_num2(nb_revealed, nb_committed):
        """ Get Upper Bound Num 2:
            Good bound for 'nb_revealed' big """
        N = nb_committed
        x = nb_revealed
        import math
        return floor((N-x)*math.log2(N/(N-x))) if x < N else 1
        
    @staticmethod
    def get_nb_leaves(nb_revealed, nb_committed, nb_experiments=1000):
        """ Get average cost to reveal the seeds """
        import random
        def run_experiment(k, N):
            logN = ceil(log2(N))
            N_ = 2**logN

            arr = [False]*N + [None]*(N_-N)
            for idx in random.sample(list(range(N)), k):
                arr[idx] = True

            count = 0
            for n in [2**i for i in range(logN, 0, -1)]:
                for i in range(0, n, 2):
                    if arr[i+1] is None:
                        arr[i//2] = arr[i]
                        continue
                    if arr[i] != arr[i+1]:
                        count += 1
                    arr[i//2] = arr[i] and arr[i+1]
                for i in range(n//2, n):
                    arr[i] = None
            count += (1 if arr[0] else 0)
            return count

        return sum([
            run_experiment(nb_revealed,nb_committed)
            for _ in range(nb_experiments)
        ]) / nb_experiments
