from framework import Search
from framework import print_title

q = 251
ws = range(-3,0+1)
nb_additional = 1/100
get_score = lambda x: 0

### Search for 128-bit security
print_title('Parameter set for 128-bit security')
print()
print('======  Hypercube variant  ======')
size, variant = Search.run(
    kappa=128,
    lda=143,
    q=q,
    n=range(220,250),
    k=range(115,140),
    w=ws,
    N=256,
    ext1=1,
    ext2=4,
    t=[3,4,5],
    nb_additional=nb_additional,
    get_score=get_score,
)
variant.print(with_sd_hardness=True,in_bytes=True)

### Search for 192-bit security
print_title('Parameter set for 192-bit security')
print()
print('======  Hypercube variant  ======')
size, variant = Search.run(
    kappa=192,
    lda=207,
    q=q,
    n=range(340,390),
    k=range(180,240),
    w=ws,
    d=2,
    N=256,
    ext1=1,
    ext2=4,
    t=[3,4,5],
    nb_additional=nb_additional,
    get_score=get_score,
)
variant.print(with_sd_hardness=True,in_bytes=True)

### Search for 256-bit security
print_title('Parameter set for 256-bit security')
print()
print('======  Hypercube variant  ======')
size, variant = Search.run(
    kappa=256,
    lda=272,
    q=q,
    n=range(460,502+1),
    k=range(250,290),
    w=ws,
    d=2,
    N=256,
    ext1=1,
    ext2=4,
    t=[3,4,5],
    nb_additional=nb_additional,
    get_score=get_score,
)
variant.print(with_sd_hardness=True,in_bytes=True)
