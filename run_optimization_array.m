function run_optimization_array(Ci)
        gpmain_rice_new(Ci)
        task_id=getenv('SLURM_ARRAY_TASK_ID');
        % Specify unique filenames for each Ci
        workspacefileName = strcat ("CO2_rice_",Ci_str,"_",task_id,".mat");
        % Save the work space 
        save(workspacefileName);
%       Save matrix of optimal enzyme rates to output file
		BestMatrix=BestMatrix'; % Transpose matrix
        BestMatrixfileName = strcat ("outputenz_",Ci_str,"_",task_id,"_rice.txt");
        d_plotfileName = strcat ("d_plot_",Ci_str,"_",task_id,"_rice.xls");
        writematrix(BestMatrix,BestMatrixfileName);
        writematrix(d_plot,d_plotfileName);
end