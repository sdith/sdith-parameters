# SD-in-the-Head: Parameter Selection

This git repository contains the scripts used to select the parameter sets for the *SD-in-the-Head signature scheme*. The latter is a digital signature scheme based on the hardness of the syndrome decoding problem for random linear codes on a finite field. SD-in-the-Head has been submitted to the NIST call for additional post-quantum signature schemes.

 * Website: [https://sdith.org/](https://sdith.org/)

## Selection Method

Extracts from the SDitH specifications:

  1. We chose to focus on syndrome decoding instances relying on the fields GF(251) and GF(256). Fields with up to 256 elements yield signature sizes close to the optimal. Moreover, for those fields an element can be stored on a byte.

  2. The remaining SD parameters (the code length m, the code dimension k and the weight w) are chosen to meet the desired security category while minimizing the signature size. For a given code length m and dimension k, the weight parameter is defined such that the number of expected solutions are below 1.01. Moreover, since the SD-in-the-Head protocol requires to have m ≤ q (to enable interpolation of a (m − 1)-degree polynomial on Fq), we use the d-split variant of the SD problem whenever necessary. In practice, we can rely on standard SD instances for Category I, and we need to rely on 2-split SD instances for Categories III and V.

  3. For the hypercube variant, we take the hypercube dimension D equal to 8 (i.e. N = 256) to achieve running times around few milliseconds while keeping short signatures. For the threshold variant, we take the maximal number N of parties allowed by the base field, namely N = q (recalling that the number of parties is at most the size of the field for this variant). The signature size is then minimized for privacy threshold ℓ = 1. However, this choice of ℓ induces an important computational overhead for commitments, thus we choose ℓ = 3 which provides a good trade-off between signature size and running times.

  4. We chose to have a common extension field for Fpoints for the two variants (hypercube and threshold) and the three security categories in order to allow a unique (optimized) implementation of the underlying extension field arithmetic. In practice, we select the degree-4 extension field Fpoints = GF(q^4) which represents a good compromise between the different settings.

  5. The parameters t and τ are chosen such that the signature size is minimal while having a forgery cost larger than λ bits, when λ is 128, 192 and 256 respectively for the categories I, III and V.

You could find more details about the selection methods in the signature specifications (available on the website).

## Usage

To run the search of the optimal parameters, you can run
```python
python3 run-search.py
```

*Remark*: the selection of the SD parameters consists in minimizing the size of the hypercube variant. The threshold variant will use the same SD parameters than the hypercube variant. For this reason, `run-search.py` only display the parameter sets for the hypercube variant. You can use the below script to get information about the threshold variant.

To just display the information (parameters, security, size, ...) of all the chosen parameter sets, we can run 
```python
python3 display-chosen.py
```

### Example

Let us take the following parameter set:
```bash
SD=(251, 235, 123, 84), nb=1+0.0086
MPC=(t=3,ext1=1,ext2=4)
Tradeoff=(N=256,tau=17) with kappa=128
```
In the above example, the displayed parameters are
 * The SD parameters `(q,m,k,w)=(251, 235, 123, 84)`, where `q` is the field size, `m` is the code length, `k` is the code dimension, and `w` is the solution weight.
 * `nb` corresponds to the averaged number of solutions for a SD instance with the given parameters.
 * The MPC parameters `(t=3,ext1=1,ext2=4)`, where `t` is the number of evaluations in the MPC protocol, `ext1` is the extension degree between the SD field and the field of the polynomial involved in the signature scheme, and `ext2` is the extension degree between the field of the polynomials and the field used for polynomial evaluations.
 * The MPCitH parameters `(N=256,tau=17)`, where `N` is the number of parties and `tau` is the number of repetitions.

## Selection Script

The selection scripts are available in the folder `framework`. Here are the description of each file:

  * `isd.py`: it contains a class `ISD` which provides several static methods to compute the cost of all the ISD algorithms for the q-ary syndrome decoding instances.
  * `sdp.py`: it contains a class `SyndromeDecoding` which represents a SD instance.
       ```python
       from framework.sdp import SyndromeDecoding
       # Create the object for the
       #   SD instance with parameters (q=256, m=230, k=126, w=79)
       sd = SyndromeDecoding(256, 230, 126, 79)
       # From then, we can get several information about the instance.
       # For example, we can get the average number of solutions
       print(sd.get_nb_solutions())
       # or its security in bits 
       print(sd.get_isd_cost())
       ```
    The class `SyndromeDecoding` also provides some useful static methods.
  * `sdith_hypercube.py`: it contains a class `HypercubeSDitH` which represents an instance of the hypercube variant of the SDitH signature.
       ```python
       from framework.sdp import SyndromeDecoding
       from framework.sdith import HypercubeSDitH
       # Create the object for the
       #   SD instance with parameters (q=256, m=230, k=126, w=79)
       sd = SyndromeDecoding(256, 230, 126, 79)
       # Create the object for the SDitH signature
       #   with parameter (t=3, ext1=1, ext2=4)
       #   for 128-bit security
       sig = HypercubeSDitH(sd, 3, 1, 4, kappa=128)
       # Set the trade-off (N=256, tau=17)
       sig.set_tradeoff(256, 17)
       # From now, we can get the signature size...
       print(sig.get_sig_size())
       # ... or its forgery cost (in bits)
       print(sig.get_signature_security())
       ```
  * `sdith_threshold.py`: it contains a class `ThresholdSDitH` which represents an instance of the threshold variant of the SDitH signature. The class `ThresholdSDitH` provides exactly the same API than `HypercubeSDitH`.
  * `search.py`: it contains a class `Search` with a unique (static) method `run`. The function `Search.run` aims to perform an exhaustive search to find the shortest signature size with the given constraints (see docstrings for details).

## Licence

This project is licensed under the terms of Apache License (version 2.0). See the [LICENSE file](LICENSE.txt) and the [copyright notice file](NOTICE).
