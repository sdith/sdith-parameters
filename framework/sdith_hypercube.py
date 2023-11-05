from math import log2, ceil, sqrt
from math import comb as binomial

class HypercubeSDitH:
    """ Represent a instance of the signature
      - sd is the syndrome decoding instance
      - t is the number of evaluations
      - ext1 is the extension degree between \F_sd and \F_poly
      - ext2 is the extension degree between \F_poly and \F_points
      - N is the number of parties
      - tau is the number of iterations
      - kappa is the security level
    """
    def __init__(self, sd, t, ext1, ext2, N=None, tau=None, kappa=128):
        self.sd = sd
        self.t = t
        self.ext1 = ext1
        self.ext2 = ext2
        self.N = N
        self.tau = tau
        self.kappa = kappa
        assert (self.sd.q**ext1) >= self.sd.n / self.sd.d
        self.p = None
        self.p = self.get_false_positive_probability()

    def set_tradeoff(self, N, tau):
        self.N = N
        self.tau = tau

    def get_parameters(self, as_tuple=False):
        if as_tuple:
            return (
                self.sd.q,
                self.sd.n,
                self.sd.k,
                self.sd.w,
                self.sd.d,
                self.t,
                self.ext1,
                self.ext2,
                self.N,
                self.tau
            )
        else:
            return {
                'q': self.sd.q,
                'n': self.sd.n,
                'k': self.sd.k,
                'w': self.sd.w,
                'd': self.sd.d,
                't': self.t,
                'ext1': self.ext1,
                'ext2': self.ext2,
                'N': self.N,
                'tau': self.tau
            }

    def get_sig_size(self):
        """ Return the signature size in bytes. More precisely, it returns (maxi, avg, std)
            where "maxi" is the maximum size, "avg" is the mean size, "std" is the standard
            deviation.
        """
        (q, _, k, w, d, t, ext1, ext2, N, tau) = self.get_parameters(as_tuple=True)

        # Components
        dig = 2*self.kappa # Digest
        salt = 2*self.kappa # Salt
        seed = self.kappa # Seed
        lN = log2(N) # Height of the generation tree
        lq = ceil(log2(q))

        # Subparts
        plaintext_size = k*lq
        poly_size = w*lq*ext1
        bn20_uni_cost = t*lq*ext1*ext2

        last_party = plaintext_size + 2*poly_size + bn20_uni_cost
        comm = 2*d*bn20_uni_cost

        # Upper bound on the signature size
        bitsize_maxi = dig + salt + tau*(
            seed*ceil(lN) + dig
            + last_party
            + comm
        )
        size_maxi = ceil(bitsize_maxi/8)

        # Mean of the signature size
        proba = (N-1)/N
        bitsize_avg = dig + salt + tau*(
            seed*ceil(lN) + dig
            + last_party*proba
            + comm
        )
        size_avg = ceil(bitsize_avg/8)

        # Standard deviation of the signature size
        #   -> This formula is valid only if n is a power of 2
        bitsize_std = last_party * sqrt(
            tau * proba * (1-proba)
        )
        size_std = ceil(bitsize_std/8)

        return size_maxi, size_avg, size_std

    @staticmethod
    def _compute_false_positive_probability(q, n, w, d, t, ext1, ext2):
        # Additional term
        delta = (q)**(ext1*ext2)
        split_n = n // d
        split_w = w // d

        pr = (split_n+split_w-1)/delta

        # Formula
        p = sum([
            binomial(t, i) * (pr)**i * (1-pr)**(t-i)
            / (delta**(t-i))
            for i in range(t+1)
        ])
        return p

    def get_false_positive_probability(self):
        if self.p is None:
            # Parameters
            (q, n, _, w, d, t, ext1, ext2, _, _) = self.get_parameters(as_tuple=True)
            self.p = self._compute_false_positive_probability(q, n, w, d, t, ext1, ext2)
        return self.p

    def get_soundness_error(self, same_randomness=False):
        (_, _, _, _, _, _, _, _, N, tau) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        if same_randomness:
            p+(1-p)*((1/N)**tau)
        else:
            (p+(1-p)/N)**tau

    @staticmethod
    def _compute_forgery_cost(p, N, tau):
        def sum_pmf(tau1, tau, p):
            return sum(
                binomial(tau, k)*(p**k)*((1-p)**(tau-k))
                for k in range(tau1, tau+1)
            )
        def inv_sum_pmf(tau1, tau, p):
            try:
                return 1/sum_pmf(tau1, tau, p)
            except ZeroDivisionError:  # Too small value
                return 2**512 # Very large value

        return log2(min(
            inv_sum_pmf(tau1, tau, p) + N**(tau-tau1)
            for tau1 in range(0, tau+1)
        ))

    def get_signature_security(self):
        """ Return the security of the signature in bits """
        (_, _, _, _, _, _, _, _, N, tau) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        return self._compute_forgery_cost(p, N, tau)

    def print(self, in_bytes=False, new_line=True, with_sd_hardness=False):
        text = []
        (q, m, k, w, d, t, ext1, ext2, n, tau) = self.get_parameters(as_tuple=True)
        p = self.get_false_positive_probability()
        nb_solutions = self.sd.get_nb_solutions()

        # SD instance
        if d == 1:
            text.append('SD=({}, {}, {}, {}), nb=1+{:.4f}'.format(q,m,k,w,nb_solutions-1))
        else:
            text.append('SD=({}, {}, {}, {}), {}-split, nb=1+{:.4f}'.format(q,m,k,w,d,nb_solutions-1))

        # MPC Protocol parameters
        text.append('MPC=(t={},ext1={},ext2={})'.format(t,ext1,ext2))

        # Trade-off
        text.append('Tradeoff=(N={},tau={}) with kappa={}'.format(n,tau,self.kappa))

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

