import sys
import numpy as np

if not len(sys.argv) == 3:
    print();
    print("\tError: Incorrect number of parameters\n");
    sys.exit();


mu = int(sys.argv[1])
sigma = int(sys.argv[2])
print(int(abs(np.random.normal(mu, sigma))))
