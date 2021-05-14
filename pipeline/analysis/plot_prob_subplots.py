import matplotlib.pyplot as plt
import tifffile as tf
import os
import matplotlib as mpl


imgs = [
    "slide_39_Ay16x1.png",
    "slide_39_Ay16x2.png",
    "slide_39_Ay16x3.png",
    "slide_16_Ay16x4.png",
    "slide_16_Ay16x5.png",
    "slide_16_Ay16x6.png",
    "slide_16_Ay16x7.png",
    "slide_16_Ay16x8.png"
    # "slide_71_Ay15x2.png",
    # "slide_71_Ay15x3.png",
    # "slide_71_Ay15x4.png",
    # "slide_71_Ay15x5.png",
    # "slide_71_Ay15x6.png",
    # "slide_39_Ay15x7.png",
    # "slide_39_Ay15x8.png"
]

# img1 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_121_39_Ay16x1_453_a0_full_maks.tiff"
# img2 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_122_39_Ay16x2_461_a0_full_maks.tiff"
# img3 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_123_39_Ay16x3_469_a0_full_maks.tiff"
# img4 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_124_16_Ay16x4_477_a0_full_maks.tiff"
# img5 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_125_16_Ay16x5_485_a0_full_maks.tiff"
# img6 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_126_16_Ay16x6_492_a0_full_maks.tiff"
# img7 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_127_16_Ay16x7_499_a0_full_maks.tiff"
# img8 = "ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_128_16_Ay16x8_506_a0_full_maks.tiff"


rows = 2
cols = 4

figs, axs = plt.subplots(rows, cols)

def get_core_name(filename: str):
    core_name = os.path.splitext(os.path.basename(filename))[0].split("_")
    if core_name[0] == "slide":
        core_name = ["ZTMA208"] + core_name
    elif core_name[0].startswith("SP"):
        core_name = ["BaselTMA"] + core_name
    
    return "_".join(core_name)
csfont = {'fontname':'Arial'}
# plt.title("Most likely \ncell type \nprobability", fontsize=7, y=-0.6, x=-3)
for r in range(rows):
    for c in range(cols):
        print(r * cols + c)
        if r * cols + c == 7: 
            break
        ax = axs[r][c]
        filename = imgs[r * cols + c]
        img = plt.imread(os.path.join("/Users/sunyunlee/Desktop/mask_imgs", imgs[r * cols + c]))
        core_name = get_core_name(filename)
        print(core_name)
        title = core_name.split("_")[-1]
        
        axs[r][c].set_title(title, fontsize=4, pad=2, **csfont)
        im = ax.imshow(img)
        ax.axis("off")


figs.delaxes(axs[1][3])
plt.subplots_adjust(wspace=0.2, hspace=-0.6, right=0.55, left=0)
# cbax = figs.add_axes([0.25, 0.23, 0.25, 0.02])
# # cbar = plt.colorbar(axes, cax = cbaxes)  


# cbar = figs.colorbar(im, ax=axs, orientation='horizontal', shrink=0.4, pad=0.05, cax=cbax, ticks=[0.4, 0.6, 0.8, 1])
# cbar.ax.tick_params(labelsize=3)
# cbar.ax.set_xlim(0.3927083311428659, 1)
# plt.clim(0.3927083311428659, 1)  # Ay16
# Ay15


ax = figs.add_axes([0.3, 0.227, 0.25, 0.02])

cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal', 
                               cmap=plt.get_cmap('viridis'), norm=mpl.colors.Normalize(0.3927083311428659, 1))
cb.ax.tick_params(labelsize=3)


plt.suptitle("Most likely cell type probability", fontsize=7, y=0.25, x=0.12, **csfont)
# plt.title("Most likely \ncell type \nprobability", fontsize=7, y=-0.6, x=-3)
# plt.savefig('plot2.png', dpi=300, bbox_inches='tight')
plt.savefig('plot_prob.pdf', dpi=1000,bbox_inches='tight')
