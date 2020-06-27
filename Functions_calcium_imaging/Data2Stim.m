function Responses =  Data2Stim(DataFold,OutFold)
 
  %%%   example : 
  %%%   DataFold = [ '~/Desktop/jacques/Data_analysis/Raw_data/',  '170821_col2/'];
  %%%   OutFold  = [ '~/Desktop/jacques/Data_analysis/Raw_data/',  '170821_col2/'];
    
    group = ''; 
    
if(~(OutFold(end)=='/')) 
    OutFold(end+1)='/'; 
end
    group = ''; 
    
if exist( [OutFold 'Data' group '.mat']  ) == 0  
    
    temp         =  find(DataFold == '/');
    DataFold1    =  DataFold(1:temp(end));
    
    load([DataFold1 'ExperimentInfo.mat'])
    load([DataFold1 savename group ' - signals.mat']);
    load([DataFold1 savename group ' - regions.mat']);

    NStimPerEp    = size(conds,1);
    Stims         = xpar.soundlist;
    Nstims        = length(Stims);
    StimDelay     = xpar.fix.StimDelay;
    TrialInterval = xpar.fix.TrialInterval;

    Data_corrected= signals -0.7*localneuropil;

    Data_corrected=  Data_corrected(7:end,:,:);
    
                                                          
 
Stimes         = ((0:NStimPerEp-1)*TrialInterval+StimDelay)/1000;
Sidx           = round(Stimes/dt);
Nbefore        = round(0.5/dt) ;   % 0.5 seconds before the stimulus
Nafter         = floor(1/dt);      % 1   s       after  


 

%%%%%%% \Delta F /F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       

%%%%%    Baseline computation ...
disp 'Compute baseline in each block';

baseline          = fn_filt(Data_corrected,50,'lk',3);
baseline          = min(baseline,[],1);
Data_corrected    = Data_corrected ./repmat(baseline,[size(Data_corrected,1) 1 1 1])-1;

%[~,  ~, drift, ~]      = spk_est(Data_corrected,'par');
%Data_corrected_v2      = (Data_corrected- drift)./drift;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Initialize and agregate data
NEp            = size(conds,2);
Ncells         = size(Data_corrected,2);
NRep           = zeros(1,Nstims);
Traces         = zeros(Ncells, Nbefore+Nafter, Nstims, floor(NEp*NStimPerEp/Nstims));
if(min(size(conds))==1)
    conds=reshape(conds,[NStimPerEp NEp]);
end

%%%%%%    DP=TransientFitting(permute(data(:,:,1:2),[2 1 3]),dt,3);

fprintf('Assigning from block %03d',0);
for i=1:NEp
    fprintf('\b\b\b%03d',i);
    Cond=conds(:,i);
    
    for j=1:NStimPerEp
       NRep(Cond(j)) = NRep(Cond(j))+1;    
                  if(NRep(Cond(j))<=floor(NEp*NStimPerEp/Nstims))            
                        Traces(:,:,Cond(j), NRep(Cond(j)))=permute(Data_corrected(Sidx(j)+(-Nbefore:Nafter-1),:,i),[2 1 3]);
                  end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    load('Definitive_corr.mat')
    Traces     = Traces(:,:,Corresp(1,:),:);
    Stims_temp = cell(size(Stims));
     for w = 1: 160
           Stims_temp{1,w}  = Stims{Corresp(1,w)};
     end
           Stims            = Stims_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

fprintf('\n');
 

[m,n,p,q]  =  size(Traces);
Calcium    = Traces - repmat( mean(Traces(:,1:10,:,:),2) ,[1 n 1 1]);                  
Calcium_av = mean(Calcium,4);



%Deconvolution and normalization

tau        =  2;  % s 
Traces_av_D= [diff(Calcium_av,1,2) zeros(m,1,p)]+[Calcium_av(:,1:end-1,:)  zeros(m,1,p)]*dt/tau;        % deconvolve

% Filter: watch out, this is not exacly what they do...

TracesDg  =gaussianfilter(Traces_av_D,1,1.5,1);

 


%%% Equation : dV/dt + V.dt/tau = firing
 
%% Saving data
if(~(OutFold(end)=='/')) 
    OutFold(end+1)='/'; 
end
%if(~exist(OutFold,'dir'))
%    mkdir(OutFold); 
%end 

  Responses             = struct;
  Responses.Traces      = Traces;
  %Responses.TracesD     = TracesD;
  %Responses.TracesDN    = TracesDN;  
  Responses.Traces_av_D = Traces_av_D;
  Responses.TracesDg    = TracesDg;     % average over repetitions, deconvoluted, filtered
  Responses.dt          = dt;
  Responses.Nbefore     = Nbefore;
  Responses.Nafter      = Nafter;
  Responses.Stims       = Stims;  
  Responses.Calcium     = Calcium;
  Responses.Calcium_av  = Calcium_av;
 
  
  save([OutFold 'Data' group '.mat'],'Traces','dt','Nbefore','Nafter','Stims','Responses')

 
 else
   
    load([OutFold 'Data' group '.mat'])

 end
  
 
  
 
  
end
  
  
  



