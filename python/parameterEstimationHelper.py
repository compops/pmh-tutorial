import numpy as np

#############################################################################
# PMH routine (LGSS)
#############################################################################
def pmh(y,initPar,par,nPart,T,xo,sm,nIter,stepSize):

    # Initalise variables
    th     = np.zeros(nIter);
    thp    = np.zeros(nIter);
    ll     = np.zeros(nIter);
    llp    = np.zeros(nIter);
    accept = np.zeros(nIter);

    # Set the initial parameter and estimate the initial log-likelihood
    th[0]  = initPar;
    ( xhat, ll[0] ) = sm(y,(th[0],par[1],par[2]),nPart,T,xo);

    #=====================================================================
    # Run main loop
    #=====================================================================
    for kk in range(1, nIter):

        # Propose a new parameter
        thp[kk] = th[kk-1] + stepSize * np.random.randn();

        # Estimate the log-likelihood
        ( xhat, llp[kk] ) = sm(y,(thp[kk],par[1],par[2]),nPart,T,xo);

        # Compute the acceptance probability
        aprob = np.min( (1.0, np.exp( llp[kk] - ll[kk-1] ) ) );

        # Generate uniform random variable in U[0,1]
        u = np.random.uniform()

        # Check if | par[0] | > 1.0, in that case set u = 1.0;
        # So that the parameter is rejected. This is due to that the model
        # is only stable when | par[0] | is smaller than 1
        if ( np.abs( thp[kk] ) > 1.0 ):
            u = 1.0;

        # Accept / reject step
        if ( u < aprob ):
            # Accept the parameter
            th[kk]     = thp[kk]
            ll[kk]     = llp[kk]
            accept[kk] = 1.0;
        else:
            # Reject the parameter
            th[kk]     = th[kk-1]
            ll[kk]     = ll[kk-1]
            accept[kk] = 0.0;

        # Write out progress
        if np.remainder(kk,100) == 0:
            print("##################################################################### ");
            print(" Iteration: " + str(kk) + " of : " + str(nIter) + " completed.")
            print("");
            print(" Current state of the Markov chain:       " + "%.4f" % th[kk] + "." )
            print(" Proposed next state of the Markov chain: " + "%.4f" % thp[kk] + "." )
            print(" Current posterior mean:                  " + "%.4f" % np.mean(th[0:kk]) + "." )
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(accept[0:kk]) + "." )
            print("##################################################################### ");
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return th;

#############################################################################
# PMH routine (SV model)
#############################################################################
def pmh_sv(y,initPar,nPart,T,sm,nIter,stepSize):

    # Initalise variables
    th     = np.zeros((nIter,3));
    thp    = np.zeros((nIter,3));
    ll     = np.zeros(nIter);
    llp    = np.zeros(nIter);
    accept = np.zeros(nIter);

    # Set the initial parameter and estimate the initial log-likelihood
    th[0,:]         = initPar;
    ( xhat, ll[0] ) = sm(y,th[0,:],nPart,T);

    #=====================================================================
    # Run main loop
    #=====================================================================
    for kk in range(1, nIter):

        # Propose a new parameter
        thp[kk,:] = th[kk-1,:] + stepSize * np.random.normal( size=3 );

        # Estimate the log-likelihood
        # Dont run if system is unstable
        if ( ( np.abs( thp[kk,0] ) < 1.0 ) & ( thp[kk,1] > 0.0 ) & ( thp[kk,2] > 0.0 ) ):
            ( xhat, llp[kk] ) = sm(y,thp[kk,:],nPart,T);

        # Compute the acceptance probability
        aprob = np.min( (1.0, np.exp( llp[kk] - ll[kk-1] ) ) );

        # Generate uniform random variable in U[0,1]
        u = np.random.uniform()

        # Check if | par[0] | > 1.0, in that case set u = 1.0;
        # Check if par[1] < 0.0, in that case set u = 1.0;
        # Check if par[2] < 0.0, in that case set u = 1.0;
        # So that the parameter is rejected. This is due to that the model
        # is only stable when | par[0] | is smaller than 1
        if ( ( np.abs( thp[kk,0] ) > 1.0 ) | ( thp[kk,1] < 0.0 ) | ( thp[kk,2] < 0.0 ) ):
            u = 1.0;

        # Accept / reject step
        if ( u < aprob ):
            # Accept the parameter
            th[kk,:]   = thp[kk,:]
            ll[kk]     = llp[kk]
            accept[kk] = 1.0;
        else:
            # Reject the parameter
            th[kk,:]   = th[kk-1,:]
            ll[kk]     = ll[kk-1]
            accept[kk] = 0.0;

        # Write out progress
        if np.remainder(kk,100) == 0:
            print("##################################################################### ");
            print(" Iteration: " + str(kk) + " of : " + str(nIter) + " completed.")
            print("");
            print(" Current state of the Markov chain:       " + "%.4f" % th[kk,0] + " " + "%.4f" % th[kk,1] + " " + "%.4f" % th[kk,2] + "." )
            print(" Proposed next state of the Markov chain: " + "%.4f" % thp[kk,0] + " " + "%.4f" % thp[kk,1] + " " + "%.4f" % thp[kk,2] + "." )
            print(" Current posterior mean:                  " + "%.4f" % np.mean(th[0:kk,0]) + " " + "%.4f" % np.mean(th[0:kk,1]) + " " + "%.4f" % np.mean(th[0:kk,2]) + "." )
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(accept[0:kk]) + "." )
            print("##################################################################### ");
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return th;
