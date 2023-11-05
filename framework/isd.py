from math import floor, log2, log
from math import comb as binomial

class ISD:
    """ Information Set Decoding Algorithm
    
        Sources:
            - Information-set decoding for linear codes over Fq,
                by Christiane Peters. https://eprint.iacr.org/2009/589.pdf
    """

    @staticmethod
    def peters_isd(n,k,q,w):
        """ Stern's adaptation of ISD over Fq, due to Peters
            It returns both complexity and optimal parameters
        """
        x = floor(k/2)

        log2q=log2(q)
        mincost=10000000
        bestp=0
        bestl=0
        max_p = min(11,floor(k/2))
        for p in range(1,max_p):
            Anum=binomial(x,p)
            Bnum=binomial(k-x,p)
            for l in range(1,floor( log(Anum)/log(q)+p*log(q-1)/log(q))+10 +1):
                ops=0.5*(n-k)**2*(n+k)+ ((0.5*k-p+1)+(Anum+Bnum)*(q-1)**p)*l+ q/(q-1.)*(w-2*p+1)*2*p*(1+(q-2)/(q-1.))*Anum*Bnum*(q-1)**(2*p)/(q**l)
                prob=Anum*Bnum*binomial(n-k-l,w-2*p)/binomial(n,w)
                cost=log2(ops)+log2(log2q)-log2(prob)
                if cost<mincost:
                    mincost=cost
                    bestp=p; bestl=l

        cost=mincost
        p=bestp
        l=bestl
        cost -= log2(q)/2
        return cost, p, l

    @staticmethod
    def lee_brickell_isd(n,k,q,w):
        """ Lee Brickell ISD over Fq
            The cost is rather approximated
                (for instance, cost of gaussian elimination is assumed to be k**3*(log2(q))**2)
        """
        p = 2
        log2_success_pr = log2(binomial(w,p)) + log2(binomial(n-w,k-p)) - log2(binomial(n,k))
        iter_cost = (n-k)**2*(n+k)*(log2(q))+binomial(k,p)*(q-1)**2*(log2(q))

        cost = log2(iter_cost)-log2_success_pr
        cost -= log2(q)/2
        return cost
