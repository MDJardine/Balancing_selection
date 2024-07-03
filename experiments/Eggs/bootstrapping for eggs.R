TLMM1.0 <- lmer(Predicted.count. ~ 1 + (1 | line) + (1 | block), data=Teggs)
TLMM1.1 <- lmer(Predicted.count. ~ allele + (1 | line) + (1 | block), data=Teggs)
TLMM1.2 <- lmer(Predicted.count. ~ allele + genotype + (1 | line) + (1 | block), data=Teggs)
TLMM1.3 <- lmer(Predicted.count. ~ allele * genotype +  (1 | line) + (1 | block), data=Teggs)

# Parametric bootstrap

# Allele

egg_0_1<-PBmodcomp(TLMM1.1, TLMM1.0)
egg_0_1
# Parametric bootstrap test; time: 45.68 sec; samples: 1000 extremes: 118;
# large : Predicted.count. ~ allele + (1 | line) + (1 | block)
# small : Predicted.count. ~ 1 + (1 | line) + (1 | block)
# stat df p.value  
# LRT    3.6342  1  0.0566 .
# PBtest 3.6342     0.1189  


# Genotype

egg_1_2<-PBmodcomp(TLMM1.2, TLMM1.1)
egg_1_2
# Parametric bootstrap test; time: 44.68 sec; samples: 1000 extremes: 43;
# large : Predicted.count. ~ allele + genotype + (1 | line) + (1 | block)
# small : Predicted.count. ~ allele + (1 | line) + (1 | block)
# stat df p.value  
# LRT    4.1691  1 0.04117 *
# PBtest 4.1691    0.04396 *


# Allele-by-Genotype interaction

egg_2_3<-PBmodcomp(TLMM1.3, TLMM1.2)
egg_2_3
# Parametric bootstrap test; time: 44.44 sec; samples: 1000 extremes: 22;
# large : Predicted.count. ~ allele * genotype + (1 | line) + (1 | block)
# small : Predicted.count. ~ allele + genotype + (1 | line) + (1 | block)
# stat df p.value  
# LRT    6.3344  1 0.01184 *
# PBtest 6.3344    0.02298 * 

