%FJ 10-2014
% Jana updated 10-2014

function initialise_CBU_Queue(userOptions)

    if userOptions.run_in_parallel || run_in_parallel_in_cluster
        try 
            matlabpool close;
        catch
            disp('Matlabpool initialising ...');
        end
        if userOptions.run_in_parallel_in_cluster
            P=cbupool;
            P.NumWorkers=userOptions.nWorkers;
            P.SubmitArguments = ['-l walltime=',num2str(userOptions.wallTime), ',mem=' ,num2str(userOptions.memReq),'gb'];      
            if isequal(userOptions.nodesReq , '^N^')
                P.ResourceTemplate = ['-l nodes=',num2str(userOptions.nodesReq)];    
            else
                P.ResourceTemplate = ['-l nodes=',num2str(userOptions.nodesReq), ':ppn=' ,num2str(userOptions.proPNode)];    
            end
            matlabpool(P);
        else
            matlabpool open; 
        end  
    end
end
