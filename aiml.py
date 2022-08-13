
import tifffile as tf
from medpy.filter import smoothing


fileName = '/home/eg/Desktop/Al-AM-samples/sample2.tif'
ofl = '/home/eg/Desktop/Al-AM-samples/output/'

image = tf.imread(fileName)
sampleName = fileName.split('/')[-1].split('.')[0]

"""Anisotropic filter
Based on trials, using the anisortropic filter in medpy, we can enhance the contrast. 

The parameters to use are: 
niter : integer, Number of iterations.

kappa : integer, Conduction coefficient, e.g. 20-100. 
   kappa controls conduction as a function of the gradient. 
   If kappa is low small intensity gradients are able to block conduction and hence diffusion across steep edges. 
   A large value reduces the influence of intensity gradients on conduction.

gamma : float, Controls the speed of diffusion. 
   Pick a value <=.25 for stability.

voxelspacing : tuple of floats or array_like. 
   The distance between adjacent pixels in all img.ndim directions

option : {1, 2, 3} Perona Malik diffusion equation No. 1 or No. 2, or Tukey’s biweight function No 3. 
   Equation 1 favours high contrast edges over low contrast ones, 
   Equation 2 favours wide regions over smaller ones.
   Equation 3 preserves sharper boundaries than previous formulations and improves the automatic stopping of the diffusion.

niter: 20
kappa: 20
gamma: 0.2
option: 3

"""

iterationVal = 20       # limit on the number of iterations
kappaVal = 20           # conduction factor
gammaVal = 20           # speed of diffusion
aisdOption = 3          # Peron-Malik vs. [Tukey] weight function

aisd_image = np.zeros_like(image)

for i in range(0, aisd_image.shape[0] + 1):
   print('Slice' + str(i))
   aisd_image[i] = smoothing.anisotropic_diffusion( img=image[i], 
                                                      voxelspacing=(1,1), 
                                                      niter=iterationVal, 
                                                      kappa=kappaVal, 
                                                      gamma=gammaVal/100, 
                                                      option=aisdOption)

saveName = sampleName + '-Iteration_' + str(iterationVal)+ '-kappa_' + str(kappaVal) + '-gamma_' + str(gammaVal) + '-Option_' + str(aisdOption) + '.tif'
tf.imwrite( ofl + saveName, aisd_image.astype('float16') )


"""Binarization with solid and void phase
"""


"""Assessment of voids as lof and other
"""


"""EXTRA code

for iteration in range( 10, 101, 10 ):
  for kappaVal in range( 20, 21 ):
      for gammaVal in range( 25, 26, 1 ):
          for option in range( 3, 4, 1 ):
              saveName = 'sample2-2-' + str(iteration)+ '-' + str(kappaVal) + '-' + str(gammaVal) + '-' + str(option) + '.tif'
              ofl = '/home/eg/Desktop/Al-AM samples/aisd/'
              print(saveName)
              aisd_image = smoothing.anisotropic_diffusion( img=image, voxelspacing=(1,1), niter=iteration, kappa=kappaVal, gamma=gammaVal/100, option=option)
              tf.imwrite( ofl + saveName, aisd_image.astype('float16') )

"""