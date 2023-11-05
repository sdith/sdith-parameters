from framework import SyndromeDecoding, HypercubeSDitH, ThresholdSDitH
from framework import print_title



print_title('Parameter set for 128-bit security')
sd_variant = SyndromeDecoding(251,242,126,87)
print()
print('======  Hypercube variant  ======')
variant = HypercubeSDitH(sd_variant, 3, 1, 4, kappa=128)
variant.set_tradeoff(256,17)
variant.print(with_sd_hardness=True, in_bytes=True)

print('======  Threshold variant  ======')
variant = ThresholdSDitH(sd_variant, 7, 1, 4, kappa=128)
variant.set_tradeoff(251,6,3)
variant.print(with_sd_hardness=True, in_bytes=True)




print_title('Parameter set for 192-bit security')
sd_variant = SyndromeDecoding(251,376,220,114,2)
print()
print('======  Hypercube variant  ======')
variant = HypercubeSDitH(sd_variant, 3, 1, 4, kappa=192)
variant.set_tradeoff(256,26)
variant.print(with_sd_hardness=True, in_bytes=True)

print('======  Threshold variant  ======')
variant = ThresholdSDitH(sd_variant, 10, 1, 4, kappa=192)
variant.set_tradeoff(251,9,3)
variant.print(with_sd_hardness=True, in_bytes=True)




print_title('Parameter set for 256-bit security')
sd_variant = SyndromeDecoding(251,494,282,156,2)
print()
print('======  Hypercube variant  ======')
variant = HypercubeSDitH(sd_variant, 4, 1, 4, kappa=256)
variant.set_tradeoff(256,34)
variant.print(with_sd_hardness=True, in_bytes=True)

print('======  Threshold variant  ======')
variant = ThresholdSDitH(sd_variant, 13, 1, 4, kappa=256)
variant.set_tradeoff(251,12,3)
variant.print(with_sd_hardness=True, in_bytes=True)
