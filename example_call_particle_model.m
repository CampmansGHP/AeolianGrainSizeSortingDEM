clear
close all



     model_run_name =   'example_generating_bed';
     disp(model_run_name);
     cd_orig = cd;
     cd(['INPUT/',model_run_name]);
            %%
            clear input_function
            [FLAGS,NUMSET,PARAMS,r,v,omega,d]=input_function;
            %%
        cd(cd_orig);
        particle_model(FLAGS,PARAMS,NUMSET,r,v,omega,d,model_run_name);
        close all;
        
        
    model_run_name = 'example_transport';
    disp(model_run_name);
    cd(['INPUT/',model_run_name]);
            %%
            clear input_function
            [FLAGS,NUMSET,PARAMS,r,v,omega,d]=input_function;
            %%
        cd(cd_orig);
        particle_model(FLAGS,PARAMS,NUMSET,r,v,omega,d,model_run_name);
        close all;
