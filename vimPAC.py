# import stuff

import numpy as np
import matplotlib.pyplot as plt

# Create variables

varA = np.zeros((5,5))
varA[3,1] = 1
varA[3,3] = 1
varA[1,1:4] = 1


# Show me what you got

print('The value of a is ',varA)

plt.figure()
plt.imshow(varA, cmap='Greys_r')
plt.xlim([0,4])
plt.ylim([0,4])
plt.show()
