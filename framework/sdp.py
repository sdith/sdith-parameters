from .isd import ISD
from math import comb as binomial
from math import log2, ceil

class SyndromeDecoding:
    """ Represent a Syndrome Decoding instance
    
      - q is the field size
      - n is the code length
      - k is the code dimension
      - w is the weight contraint
      - d is the split factor (d=1 for standard SD instance)
    """

    def __init__(self, q, n, k, w, d=1):
        self.q = q
        self.n = n
        self.k = k
        self.w = w
        self.d = d

        # cache
        self._cost_peters_isd = None
        self._cost_lee_brickell_isd = None
        assert (n % d == 0) and (w % d == 0)

    def get_security_loss_from_split(self):
        from math import log2
        from math import comb as binomial
        split_n = self.n // self.d
        split_w = self.w // self.d
        return log2(binomial(self.n,self.w)) - log2(binomial(split_n,split_w)**self.d)

    @staticmethod
    def compute_max_weigth_for_target(q, n, k, ratio=1/100):
        """ This function returns the maximal weight w
            such that there is at most 1+ratio solutions
            in average to the syndrome decoding problem.

            The averaged number of solutions is
                1 + (sum_{i=0}^w binom(n,i)*(q-1)**i)/(q**(n-k)).
            So, it want the maximal w such that
                (sum_{i=0}^w binom(n,i)*(q-1)**i)/(q**(n-k)) <= ratio.
            By writing
                left_term := (sum_{i=0}^w binom(n,i)*(q-1)**i)/(q**(n-k))
                right_term := q**(n-k),
            we get that we search the maximal w such that
                left_term / right_term <= ratio,
            or equivently,
                left_term * (1/ratio) <= right_term.
        """
        right_term = q**(n-k) # target

        d = 0
        left_term = 1
        # Python does not support working on very large numbers
        #   represented as floats, so we added `ceil` to
        #   work over integers.
        while left_term*ceil(1/ratio) <= right_term:
            # Loop Invariant:
            #   left_term = sum_{i=0}^{d} binom(n,i)*(q-1)**i
            d+=1
            left_term += binomial(n,d)*(q-1)**d

        # The current `d` is the first `d`
        #   such that left_term / right_term > ratio,
        #   so we must take the previous one.
        d = d-1
        return d

    def get_max_weight_for_target(self, ratio=1/100):
        q = self.q
        n = self.n
        k = self.k
        return self.compute_max_weigth_for_target(q, n, k, ratio)

    @staticmethod
    def compute_nb_solutions(q, n, k, w):
        """ The averaged number of solutions of a SD problem is
                1 + (sum_{i=0}^w binom(n,i)*(q-1)**i)/(q**(n-k)).
        """
        right_term = q**(n-k)
        left_term = sum(
            binomial(n,d)*(q-1)**d
            for d in range(w+1)
        )
        return 1+left_term/right_term

    def get_nb_solutions(self):
        q = self.q
        n = self.n
        k = self.k
        w = self.w
        return self.compute_nb_solutions(q, n, k, w)
    
    def get_cost_peters_isd(self, with_parameters=False):
        if self._cost_peters_isd is None:
            self._cost_peters_isd = ISD.peters_isd(self.n,self.k,self.q,self.w)
        cost, p, ell = self._cost_peters_isd
        if with_parameters:
            return cost, (p, ell)
        return cost - self.get_security_loss_from_split()

    def get_cost_lee_brickell_isd(self):
        if self._cost_lee_brickell_isd is None:
            self._cost_lee_brickell_isd = ISD.lee_brickell_isd(self.n,self.k,self.q,self.w)
        cost = self._cost_lee_brickell_isd
        return cost - self.get_security_loss_from_split()

    def get_isd_cost(self):
        cost_peters_isd = self.get_cost_peters_isd()
        cost_lee_brickell_isd = self.get_cost_lee_brickell_isd()
        return min(cost_peters_isd, cost_lee_brickell_isd)
