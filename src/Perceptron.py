import numpy as np 

class Perceptron:
    
    def __init__(self):
        pass
    
    @staticmethod
    def predict(weights, instance):
        
        return 1 if weights.dot(instance) > 0 else 0
        
    
    def fit(self, instances, labels, epochs= 10, alpha= .01, verbose= True):
        
        n = instances.shape[1]
        
        weights = np.zeros([1, n])
    
        for epoch in range(epochs):
            error = 0
            
            for instance, y in zip(instances, labels):
            
                y_ = self.predict(weights, instance)
            
            # update weights
            
                weights += alpha * (y - y_) * instance
            
                if y != y_:
                    error += 1
                
            error_rate = error / len(instances)
            
            if verbose:
                print("[%d] \t Error rate: %f" % (epoch, error_rate))
        
        return weights
        
        
    

if __name__ == "__main__":

    import seaborn as sns
    import matplotlib.pyplot as plt 

    x1 = np.random.normal(loc=0, scale=1, size=100)
    x2 = np.random.normal(loc=0, scale=1, size=100)

    y  = np.where(x1+x2 > 0, 1, 0)

    print('First see the scatter plot')
    sns.scatterplot(x1, x2, hue= y)
    plt.show()
    print('Start training process ...') 

    X = np.c_[x1, x2]

    tron = Perceptron()

    weight = tron.fit(X, y)

    print(f"The weight is {weight}")

    print("Scatter plot:")

    sns.scatterplot(x1, x2, hue=y)
    x_plot = np.linspace(-2, 2, 100)
    y_plot = -(weight[0][0] / weight[0][1]) * x_plot
    plt.plot(x_plot, y_plot, color='r')
    plt.show()
    
    
    
    
    
    