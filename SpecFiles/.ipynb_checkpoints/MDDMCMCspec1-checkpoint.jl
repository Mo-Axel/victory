nsim       = 5000
nburn      = 5000

l1_grid    = exp.(range(-7, stop=6, length=10))
l1_none    = ones(size(l1_grid))
l2_grid    = exp.(range(-7, stop=6, length=10))
l2_none    = ones(size(l2_grid))
l3_grid    = exp.(range(-7, stop=6, length=10))
l3_none    = ones(size(l3_grid))

# l1_grid    = exp.(range(-7, stop=6, length=15))
# l1_none    = ones(size(l1_grid))
# l2_grid    = exp.(range(-7, stop=6, length=15))
# l2_none    = ones(size(l2_grid))
# l3_grid    = exp.(range(-7, stop=6, length=15))
# l3_none    = ones(size(l3_grid))

hyper_grid = [ kron(l1_grid,l2_none,l3_none) kron(l1_none,l2_grid,l3_none) kron(l1_none,l2_none,l3_grid)]
hyper_n    = size(hyper_grid)[1]
MDD_term2_vec  = zeros(hyper_n);

tau     = 0.05 # for Geweke's MDD approximation
seedoffset = 0
