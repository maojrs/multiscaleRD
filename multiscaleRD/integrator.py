class integrator:   

    def LangevinSolution(setting, particle, fluid):
        T = setting.intervallend
        g = particle.potentialenergy
        gamma = fluid.dampingcoefficient
        N = setting.Brownian_paths
        D = fluid.constant
        dt=T/N
        
        np.random.seed( N )
        dW = np.sqrt(T/N)*np.random.randn(1,N)
    
        t = np.linspace(dt, T, N)
        dW2 = np.concatenate(dW)
        W = np.cumsum(dW2)
        
        R = setting.R
        Dt = R*dt
        L = int(N/R)
        
        Xem = numpy.zeros(L)
        Xem2 = numpy.zeros(L)
        Xzero = particle.position
        Xem[0] = Xzero
        Xem2[0] = Xzero
        j = 1
        
       
        
        while j < L:
         
            #Xem[j] = Xem[j-1]+Dt/gamma*(V(Xem[j-1])-V(Xem[j-2]))/(Xem[j-1]-Xem[j-2])+np.sqrt(2*D)*np.sum(dW2[R*(j-1)+1:R*j])
            Xem2[j]= Xem2[j-1]+Dt/gamma*g(Xem2[j-1])+np.sqrt(2*D)*np.sum(dW2[R*(j-1)+1:R*j])
            j = j+1
       
            
        
        #plt.plot(np.linspace(0,T, int(T/(Dt))), Xem, 'r')
        plt.plot(np.linspace(0,T, int(T/(Dt))), Xem2, 'b')
    
    def LangevinSolution2D(setting, particle, fluid):
        
        
        
    
        
        T = setting.intervallend
        g = particle.potentialenergy
        gamma = fluid.dampingcoefficient
        N = setting.Brownian_paths
        D = fluid.constant
        R = setting.R
        Xzero = particle.position
                    
        dt=T/N
        np.random.seed( N )
        dW = np.sqrt(T/N)*np.random.randn(2,N)# 2 dimensions
               
                #t = np.linspace(dt, T, N)
                
        W = numpy.zeros(shape=(2,N))
        W[0,:] = np.cumsum(dW[0,:])
                
        W[1,:] = np.cumsum(dW[1,:])
                
        Dt = R*dt
        L = int(N/R)
                
                
        Xem2 = numpy.zeros(shape=(2,L))
                
        Xem2[:,0] = Xzero
                
               
        j = 1
                
               
                
        while j in range(L):
                 
                    #Xem[0,j] = Xem[0,j-1]+Dt/gamma*(V(Xem[0,j-1])-V(Xem[0,j-2]))/(Xem[0,j-1]-Xem[0,j-2])+np.sqrt(2*D)*np.sum(dW2[R*(j-1)+1:R*j])
                    #Xem[1,j] = Xem[1,j-1]+Dt/gamma*(V(Xem[1,j-1])-V(Xem[1,j-2]))/(Xem[1,j-1]-Xem[1,j-2])+np.sqrt(2*D)*np.sum(dW2[R*(j-1)+1:R*j])
            force = Dt/gamma*g(Xem2[:,j-1])
            Xem2[0,j]= Xem2[0,j-1]-force[0]+np.sqrt(2*D)*np.sum(dW[0,:][R*(j-1)+1:R*j])
                
            Xem2[1,j]= Xem2[1,j-1]-force[1]+np.sqrt(2*D)*np.sum(dW[1,:][R*(j-1)+1:R*j])
            j = j+1 
                  
                
        plt.plot(Xem2[0,:],Xem2[1,:], 'r')
        fig = plt.figure()
                
            
        ax = plt.axes(projection='3d')
        ax.plot(np.linspace(0,T, int(T/(Dt))),Xem2[0,:],Xem2[1,:],'blue' )   


    
        
    def Movement(setting, particle):
        
        T = setting.intervallend
        N= setting.Brownian_paths
        M = setting.M
        f = particle.f
        
        np.random.seed( 2 )
        dt = T/N
        t = np.linspace(dt, T, N)
        dW = np.sqrt(dt)*np.random.randn(M,N)
           
        W = numpy.zeros(shape=(M,N))
        F = numpy.zeros(shape=(M,N))
        for j  in range(M):
            F[j,:] = f(np.cumsum(dW[j,:]))
            plt.plot(t, F[j,:])
            plt.xlabel('t')
            plt.ylabel('f(W)')
        
                

    
    def EulerMaruyumaSolution(setting, particle):
        T = setting.intervallend
        g = particle.potentialenergy
        f = particle.f
        
        N = setting.Brownian_paths
        
        R = setting.R
        Xzero = particle.position
        
        np.random.seed( N )
        dt=T/N
        dW = np.sqrt(T/N)*np.random.randn(1,N)
        dW2 = np.concatenate(dW)
        t = np.linspace(dt, T, N)
        W = np.cumsum(dW2)
        
        Dt = R*dt
        L = int(N/R)
        
        Xem = numpy.zeros(L)
        Xem[0] = Xzero
        j = 1
        
       
       
        while j < L:
         
            
          
            
               
                
            Xem[j] = Xem[j-1]+Dt*f(Xem[j-1])+g(Xem[j-1])*np.sum(dW2[R*(j-1)+1:R*j])
            j = j+1
        
        r = np.exp((-t)+2*W) #Testfunktion
        
        
        
        plt.plot(t, r) # richtiger Wert
        plt.plot(np.linspace(0,T, int(T/(Dt))), Xem, 'r')