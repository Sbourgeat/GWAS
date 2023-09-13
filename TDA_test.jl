using Ripserer
using Images
using Plots
using PersistenceDiagrams

fly_brain = load(
    "/Users/skumar/TDA/Cubical_complex/only_brain/32_female_all.tif"
)

plot(fly_brain[:,:,150])

result = ripserer(Cubical(-fly_brain), dim_max=2)

plot(result)


