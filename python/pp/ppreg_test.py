import ppreg
import matplotlib.pyplot as plt
import numpy as np
       
if __name__ == "__main__":
    n = 200

    x1 = np.random.uniform(-1,1,size=n)
    x2 = np.random.uniform(-1,1,size=n)
    error = np.random.normal(scale=0.04,size=n)
    y = x1*x2 + error

    x = np.empty((n,2))
    x[:,0] = x1
    x[:,1] = x2

    model = ppreg.ProjectionPursuitRegression(3,1)
    model.fit(x,y)

    prediction = model.predict(x)

    print np.corrcoef(y,prediction)

    proj = np.dot(x,model.alpha[0])

    plt.scatter(proj,y)

    plt.show()

