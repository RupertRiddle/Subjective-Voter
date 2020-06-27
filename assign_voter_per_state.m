
opts = detectImportOptions('2012_test');

state_data = readtable('2012_test',opts);
tv = 99; % 1-N: sets the number of voters
Ny = 1; % N+1 : Sets the number of election cycles not including the first election
win = 2; % 1-2: sets the initial party in power for all agents
good_econ = 0.8; %how much of the economy is good
bad_econ = 0.2; % how much of the economy is bad
alpha = 2.33; % sets the shape of the distribution for preferences
beta = 1;        
swing_precision_lb = 1;
swing_precision_ub = 2;


Partisanship = 1.034483; %volatility of mappings for republicans(can be used to simulate extremeness)
swing_polarisation_democrat_bias = 1;
swing_polarisation_republican_bias = 1.11;




% Markov Chain Monte Carlo generates distribution of prior beliefs and
% preferences
victory_margin = [];
pdf = @(x)gampdf(x,alpha,beta); % Target distribution
proppdf = @(x,y)gampdf(x,floor(alpha),floor(alpha)/alpha);
proprnd = @(x)sum(...
              exprnd(floor(alpha)/alpha,floor(alpha),1));
mdp_array = [];
nsamples = tv;
smpl = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
len = length(smpl);


winhist=zeros(Ny+1,1);

winhist(1,1)=win;
total_states = 51;
    for i =1:total_states
        row = state_data(i,:);
        state = table2array(row(:,1));
        n_rep = table2array(row(:,3));
        n_dem = table2array(row(:,2));
        n_und = table2array(row(:,4));
        n__dem = floor(tv*(floor(n_dem)/100));
        n__rep = floor(tv*(floor(n_rep)/100));
        n__und = tv - (n__dem + n__rep);
        
        for nv = 1:tv
            if nv <= n__rep
                for i = 1:3
                    econ_index = randi([1 len],1,1); % generates random number for indexing smpl probability distribution
                    social_index = randi([1 len],1,1);
                    elec_index= randi([1 len],1,1);
                    econ_pref = smpl(econ_index);
                    social_pref = smpl(social_index);
                    elec_pref = smpl(elec_index);
                    precision_index = randi([1 len],1,1);

                end
                T=2; %Econonmy, Social, Election, Satisfaction 
        %%%republican:
        %%%if vote = rep,  (voting impact beliefs)
        %%%-> economy = good, social=right, election=rep
        %%%if power = rep, (government impact beliefs)
        %%%-> economy = good, social=right, election=rep
        %%%if satisfaction = good, (satisfaction inference)
        %%%-> economy = good, social=right, election=rep


                A{1} = zeros(2,3,2,2) +0.5; %economy outcomes: good/bad, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{1}(2,:,1,:) = A{1}(2,:,1,:)*(Partisanship*3.0); %%dem power --> bad economy
                A{1}(1,:,2,:) = A{1}(1,:,2,:)*(Partisanship*3.0); %rep power -> Good Economy
                A{1}(1,:,:,1) = A{1}(1,:,:,1)*(Partisanship*3.0); %%satisfaction --> good economy
                A{1}(2,:,:,2) = A{1}(2,:,:,2)*(Partisanship*3.0); %%dissatisfaction --> bad economy

                A{2} = zeros(2,3,2,2) +0.5; %social outcomes: left/right, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{2}(2,:,2,:) = A{1}(2,:,2,:)*(Partisanship*3.0); %%rep power --> right-oriented
                A{2}(1,:,1,:) = A{1}(1,:,1,:)*(Partisanship*3.0); %%dem power --> left-oriented
                A{2}(2,:,:,1) = A{1}(2,:,:,1)*(Partisanship*3.0); %%satisfaction --> right-oriented
                A{2}(1,:,:,2) = A{1}(1,:,:,2)*(Partisanship*3.0); %%dissatisfaction --> left-oriented


                A{3} = zeros(2,3,2,2) +0.5; %election outcomes: dem/rep, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{3}(2,2,:,:) = A{3}(2,2,:,:)*(Partisanship*1.5); %%% vote --> election outcome
                A{3}(1,1,:,:) = A{3}(1,1,:,:)*(Partisanship*1.5);

                A{3}(2,:,2,:) = A{3}(2,:,2,:)*100.0; %%%power state ->> election outcome
                A{3}(1,:,1,:) = A{3}(1,:,1,:)*100.0; 
                beta_A = 2;
                for i = 1:2; for j = 1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=A{1}(i,j,k,l)^(beta_A^-1);
                    An{2}(i,j,k,l)=A{2}(i,j,k,l)^(beta_A^-1);
                    An{3}(i,j,k,l)=A{3}(i,j,k,l)^(beta_A^-1);
                end;end;end;end

                normAn = zeros(3,2,2);
                normAn2 = zeros(3,2,2);
                normAn3 = zeros(3,2,2);

                for j=1:3; for k=1:2; for l = 1:2
                    normAn(j,k,l) = sum(An{1}(:,j,k,l));
                    normAn2(j,k,l) = sum(An{2}(:,j,k,l));
                    normAn3(j,k,l) = sum(An{3}(:,j,k,l));
                end;end;end

                for i =1:2; for j=1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=An{1}(i,j,k,l)*normAn(j,k,l)^-1;
                    An{2}(i,j,k,l)=An{2}(i,j,k,l)*normAn2(j,k,l)^-1;
                    An{3}(i,j,k,l)=An{3}(i,j,k,l)*normAn3(j,k,l)^-1;
                end;end;end;end


                %transition
                B{1} = zeros(3,3,3); %state transitions for vote (dem/rep/none) to vote (dem/rep/none)
                for i = 1:3
                     B{1}(i,:,i) = 1;
                end


                B{2} = zeros(2,2,3)+0.5; %%%two transitions, one if vote = dem, one if vote = rep
                B{2}(1,:,1) = B{2}(1,:,1)*(Partisanship*1.5); %%vote dem --> transition to dem power more likely
                B{2}(2,:,2) = B{2}(2,:,2)*(Partisanship*1.5); %%vote rep --> transition to rep power more likely

                B{3} = zeros(2,2,3)+0.5;
                B{3}(2,:,1) = B{3}(2,:,1)*(Partisanship*1.5); %%vote dem --> transition to dissatisfaction more likely
                B{3}(1,:,2) = B{3}(1,:,2)*(Partisanship*1.5); %%vote rep --> transition to satisfaction more likely


                % allowable actions
                U = zeros(1,3,3); %policy matrix (1 time step)
                U(1,:,1)=1:3; %%vote dem/rep/none
                U(1,:,2)=1; %%
                U(1,:,3)=1;

                % preferences over outcomes 
                C{1}=zeros(2,T); %economy outcomes
                C{1}(1,:) = econ_pref; %prefers good outcome
                C{1}(2,:) = -econ_pref; %fears bad economy

                C{2}=zeros(2,T); %social outcomes
                C{2}(1,:) = -social_pref; %prefers right
                C{2}(2,:) = social_pref; %dislikes left

                C{3}=zeros(2,T); %election outcomes
                C{3}(1,:) = -elec_pref; %prefers republican
                C{3}(2,:) = elec_pref; %dislikes democrat

                D{1} = zeros(3,1); %initial vote states (dem/rep/none)
                D{1}(2,1) = 1.0;

                D{2} = zeros(2,1); %power states (dem/rep)
                D{2} = [0.5,0.5]';

                D{3} = zeros(2,1); %satisfaction states (+,-)
                D{3} = [0.5,0.5]';

                E = zeros(3,1);   
                E= [0.1,0.7,0.2]';
                label.factor{1}   = 'vote'; label.name{1}    = {'dem','rep','none'};
                label.factor{2}   = 'power';  label.name{2}    = {'dem','rep'};
                label.factor{3}   = 'satisfaction';     label.name{3}    = {'pos','neg'};
                label.action{1} = {'dem','rep','none'};
                label.modality{1} ='economy'; label.outcome{1}={'good','bad'};
                label.modality{2} ='social'; label.outcome{2}={'left','right'};
                label.modality{3} ='election'; label.outcome{3}={'dem','rep'};
                %mdp(nv,1).o(1,:)=[econ];
                mdp(nv,1).O{1}(1,:)= good_econ;
                mdp(nv,1).O{1}(2,:)= bad_econ;
                mdp(nv,1).o(2,:)=[win]';
                mdp(nv,1).o(3,:)=[win]';

                mdp(nv,1).label = label;
                mdp(nv,1).T =T;
                mdp(nv,1).a =An;
                mdp(nv,1).B =B;
                mdp(nv,1).C =C;
                mdp(nv,1).D =D;
                mdp(nv,1).E =E;
                %mdp(nv,1).beta = 1.0;
                mdp(nv,1).U =U;
                %spm_MDP_check(mdp)
                MDP(nv,1) = spm_MDP_VB_X(mdp(nv,1));
                %spm_MDP_VB_trial(MDP(nv,1)); taken out so it doesnt generate a
                %graph
                %sat = MDP(nv,1).X{3}(:,2)

                mdp_array = [mdp_array, nv];
            elseif (nv <= (n__dem + n__rep) && (nv > n__rep))

                len = length(smpl);
                for i = 1:3
                    econ_index = randi([1 len],1,1); % generates random number for indexing smpl probability distribution
                    social_index = randi([1 len],1,1);
                    elec_index= randi([1 len],1,1);
                    econ_pref = smpl(econ_index);
                    social_pref = smpl(social_index);
                    elec_pref = smpl(elec_index);
                    precision_index = randi([1 len],1,1);

                end
                T = 2;

                A{1} = zeros(2,3,2,2) +0.5; %economy outcomes: good/bad, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{1}(1,:,1,:) = A{1}(1,:,1,:)*(Partisanship *3.0); %%dem power --> good economy
                A{1}(2,:,2,:) = A{1}(2,:,2,:)*(Partisanship *3.0); %%rep power --> bad economy
                % Party in power is penalised for bad economy. A{1}


                A{1}(1,:,:,1) = A{1}(1,:,:,1)*(Partisanship *3.0); %%satisfaction --> good economy
                A{1}(2,:,:,2) = A{1}(2,:,:,2)*(Partisanship *3.0); %%dissatisfaction --> bad economy

                A{2} = zeros(2,3,2,2) +0.5; %social outcomes: left/right, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...

                A{2}(2,:,2,:) = A{1}(2,:,2,:)*(Partisanship *3.0); %%rep power --> right-oriented
                A{2}(1,:,1,:) = A{1}(1,:,1,:)*(Partisanship *3.0); %%dem power --> left-oriented
                A{2}(2,:,:,2) = A{1}(2,:,:,2)*(Partisanship *3.0); %%satisfaction --> right-oriented
                A{2}(1,:,:,1) = A{1}(1,:,:,1)*(Partisanship *3.0); %%dissatisfaction --> left-oriented


                A{3} = zeros(2,3,2,2) +0.5; %election outcomes: dem/rep, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{3}(2,2,:,:) = A{3}(2,2,:,:)*(Partisanship *1.5); %%% vote --> election outcome
                A{3}(1,1,:,:) = A{3}(1,1,:,:)*(Partisanship *1.5);

                A{3}(2,:,2,:) = A{3}(2,:,2,:)*100.0; %%%power state ->> election outcome
                A{3}(1,:,1,:) = A{3}(1,:,1,:)*100.0; 
                beta_A = 2;
                for i = 1:2; for j = 1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=A{1}(i,j,k,l)^(beta_A^-1);
                    An{2}(i,j,k,l)=A{2}(i,j,k,l)^(beta_A^-1);
                    An{3}(i,j,k,l)=A{3}(i,j,k,l)^(beta_A^-1);
                end;end;end;end

                normAn = zeros(3,2,2);
                normAn2 = zeros(3,2,2);
                normAn3 = zeros(3,2,2);

                for j=1:3; for k=1:2; for l = 1:2
                    normAn(j,k,l) = sum(An{1}(:,j,k,l));
                    normAn2(j,k,l) = sum(An{2}(:,j,k,l));
                    normAn3(j,k,l) = sum(An{3}(:,j,k,l));
                end;end;end

                for i =1:2; for j=1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=An{1}(i,j,k,l)*normAn(j,k,l)^-1;
                    An{2}(i,j,k,l)=An{2}(i,j,k,l)*normAn2(j,k,l)^-1;
                    An{3}(i,j,k,l)=An{3}(i,j,k,l)*normAn3(j,k,l)^-1;
                end;end;end;end



                %transition
                B{1} = zeros(3,3,3); %state transitions for vote (dem/rep/none) to vote (dem/rep/none)
                for i = 1:3
                   B{1}(i,:,i) = 1;
                end
                B{2} = zeros(2,2,3)+0.5; %%%two transitions, one if vote = dem, one if vote = rep
                B{2}(1,:,1) = B{2}(1,:,1)*(Partisanship *1.5); %%vote dem --> transition to dem power more likely
                B{2}(2,:,2) = B{2}(2,:,2)*(Partisanship *1.5); %%vote rep --> transition to rep power more likely

                B{3} = zeros(2,2,3)+0.5;
                B{3}(2,:,2) = B{3}(2,:,2)*(Partisanship *1.5); %%vote rep --> transition to dissatisfaction more likely
                B{3}(1,:,1) = B{3}(1,:,1)*(Partisanship *1.5); %%vote dem --> transition to satisfaction more likely


                % allowable actions
                U = zeros(1,3,3); %policy matrix (1 time step)
                U(1,:,1)=1:3; %%vote dem/rep/none
                U(1,:,2)=1; %%
                U(1,:,3)=1;

                % preferences over outcomes 
                C{1}=zeros(2,T); %economy outcomes
                C{1}(1,:) = econ_pref; %prefers good outcome
                C{1}(2,:) = -econ_pref; %fears bad economy

                C{2}=zeros(2,T); %social outcomes
                C{2}(1,:) = social_pref; %prefers left
                C{2}(2,:) = -social_pref; %dislikes right

                C{3}=zeros(2,T); %election outcomes
                C{3}(1,:) = elec_pref; %prefers dem
                C{3}(2,:) = -elec_pref; %dislikes rep

                D{1} = zeros(3,1); %initial vote states (dem/rep/none)
                D{1}(2,1) = 1.0;

                D{2} = zeros(2,1); %power states (dem/rep)
                D{2} = [0.5,0.5]';

                D{3} = zeros(2,1); %satisfaction states (+,-)
                D{3} = [0.5,0.5]';

                E = zeros(3,1);   
                E= [0.7,0.1,0.2]';
                label.factor{1}   = 'vote'; label.name{1}    = {'dem','rep','none'};
                label.factor{2}   = 'power';  label.name{2}    = {'dem','rep'};
                label.factor{3}   = 'satisfaction';     label.name{3}    = {'pos','neg'};
                label.action{1} = {'dem','rep','none'};
                label.modality{1} ='economy'; label.outcome{1}={'good','bad'};
                label.modality{2} ='social'; label.outcome{2}={'left','right'};
                label.modality{3} ='election'; label.outcome{3}={'dem','rep'};

                mdp(nv,1).O{1}(1,:)= good_econ;
                mdp(nv,1).O{1}(2,:)= bad_econ;
                mdp(nv,1).o(2,:)=[win]'; %how is the social policy
                mdp(nv,1).o(3,:)=[win]'; %what is the election result

                mdp(nv,1).label = label;
                mdp(nv,1).T =T;
                mdp(nv,1).a =An;
                mdp(nv,1).B =B;
                mdp(nv,1).C =C;
                mdp(nv,1).D =D;
                mdp(nv,1).E =E;
                %mdp(1,1).beta = 1.0; (low precision:beta=1.5, high precision:beta=0.5)
                mdp(nv,1).U =U;
                %spm_MDP_check(mdp)
                MDP(nv,1) = spm_MDP_VB_X(mdp(nv,1));
                mdp_array = [mdp_array, nv];


            else
                %Insert Swing Voterfor i = 1:3
                for i= 1:2
                    dem_index = randi([1 len],1,1); % generates random number for indexing smpl probability distribution
                    rep_index= randi([1 len],1,1);
                    dem_belief = smpl(dem_index);
                    rep_belief = smpl(rep_index);
                    precision_index = randi([1 len],1,1);
                    
                end
                T = 2;
                A{1} = zeros(2,3,2,2) +0.5; %economy outcomes: good/bad, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{1}(2,:,1,2) = A{1}(2,:,1,2) *(swing_polarisation_democrat_bias*3.0);
                A{1}(2,:,2,2) = A{1}(2,:,2,2) *(swing_polarisation_republican_bias*3.0);
                % Party in power is rewarded for Good economy. 
                A{1}(2,:,1,1) = A{1}(2,:,1,1) *(swing_polarisation_democrat_bias*3.0);
                A{1}(2,:,2,1) = A{1}(2,:,2,1) *(swing_polarisation_republican_bias*3.0);

                % Randomised belief 
                A{1}(1,:,1,:) = A{1}(1,:,1,:)*(swing_polarisation_democrat_bias*dem_belief); %%dem power --> Good economy
                A{1}(1,:,2,:) = A{1}(1,:,2,:)*(swing_polarisation_republican_bias*rep_belief); %rep power -> Good Economy

                A{1}(1,:,:,1) = A{1}(1,:,:,1)*3.0; %%satisfaction --> good economy
                A{1}(2,:,:,2) = A{1}(2,:,:,2)*3.0; %%dissatisfaction --> bad economy

                A{2} = zeros(2,3,2,2) +0.5; %social outcomes: left/right, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...

                A{2}(2,:,2,:) = A{1}(2,:,2,:)*3.0; %%rep power --> right-oriented
                A{2}(1,:,1,:) = A{1}(1,:,1,:)*3.0; %%dem power --> left-oriented

                A{3} = zeros(2,3,2,2) +0.5; %election outcomes: dem/rep, states: vote (dem/rep/none), power (dem/rep), sat(+,-): Likely observation of a...
                A{3}(2,2,:,:) = A{3}(2,2,:,:)*1.5; %%% vote --> election outcome
                A{3}(1,1,:,:) = A{3}(1,1,:,:)*1.5;

                A{3}(2,:,2,:) = A{3}(2,:,2,:)*100.0; %%%power state ->> election outcome
                A{3}(1,:,1,:) = A{3}(1,:,1,:)*100.0; 

                r = (swing_precision_ub-swing_precision_lb).*rand(tv,1) + swing_precision_lb;
                beta_A = r(precision_index);
                for i = 1:2; for j = 1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=A{1}(i,j,k,l)^(beta_A^-1);
                    An{2}(i,j,k,l)=A{2}(i,j,k,l)^(beta_A^-1);
                    An{3}(i,j,k,l)=A{3}(i,j,k,l)^(beta_A^-1);
                end;end;end;end

                normAn = zeros(3,2,2);
                normAn2 = zeros(3,2,2);
                normAn3 = zeros(3,2,2);

                for j=1:3; for k=1:2; for l = 1:2
                    normAn(j,k,l) = sum(An{1}(:,j,k,l));
                    normAn2(j,k,l) = sum(An{2}(:,j,k,l));
                    normAn3(j,k,l) = sum(An{3}(:,j,k,l));
                end;end;end

                for i =1:2; for j=1:3; for k = 1:2; for l = 1:2
                    An{1}(i,j,k,l)=An{1}(i,j,k,l)*normAn(j,k,l)^-1;
                    An{2}(i,j,k,l)=An{2}(i,j,k,l)*normAn2(j,k,l)^-1;
                    An{3}(i,j,k,l)=An{3}(i,j,k,l)*normAn3(j,k,l)^-1;
                end;end;end;end



                %transition
                B{1} = zeros(3,3,3); %state transitions for vote (dem/rep/none) to vote (dem/rep/none)
                for i = 1:3
                   B{1}(i,:,i) = 1;
                end


                B{2} = zeros(2,2,3)+0.5; %%%two transitions, one if vote = dem, one if vote = rep
                B{2}(1,:,1) = B{2}(1,:,1)*1.5; %%vote dem --> transition to dem power more likely
                B{2}(2,:,2) = B{2}(2,:,2)*1.5; %%vote rep --> transition to rep power more likely

                B{3} = zeros(2,2,3)+0.5;



                % allowable actions
                U = zeros(1,3,3); %policy matrix (1 time step)
                U(1,:,1)=1:3; %%vote dem/rep/none
                U(1,:,2)=1; %%
                U(1,:,3)=1;

                % preferences over outcomes 
                C{1}=zeros(2,T); %economy outcomes
                C{1}(1,:) = 5; %prefers good outcome
                C{1}(2,:) = -5; %fears bad economy

                C{2}=zeros(2,T); %social outcomes
                C{2}(1,:) = 0; %prefers right
                C{2}(2,:) = 0; %dislikes left

                C{3}=zeros(2,T); %election outcomes
                C{3}(1,:) = 0; %prefers republican
                C{3}(2,:) = 0; %dislikes democrat

                D{1} = zeros(3,1); %initial vote states (dem/rep/none)
                D{1}(2,1) = 3.0;


                D{2} = zeros(2,1); %power states (dem/rep)
                D{2} = [0.5,0.5]';

                D{3} = zeros(2,1); %satisfaction states (+,-)
                D{3} = [0.5,0.5]';

                E = zeros(3,1);   
                E= [0.45,0.45,0.05]';
                mdp(nv,1).O{1}(1,:)= good_econ;
                mdp(nv,1).O{1}(2,:)= bad_econ;
                mdp(nv,1).o(2,:)=[win]'; %how is the social policy
                mdp(nv,1).o(3,:)=[win]'; %what is the election result
                label.factor{1}   = 'vote'; label.name{1}    = {'dem','rep','none'};
                label.factor{2}   = 'power';  label.name{2}    = {'dem','rep'};
                label.factor{3}   = 'satisfaction';     label.name{3}    = {'pos','neg'};
                label.action{1} = {'dem','rep','none'};
                label.modality{1} ='economy'; label.outcome{1}={'good','bad'};
                label.modality{2} ='social'; label.outcome{2}={'left','right'};
                label.modality{3} ='election'; label.outcome{3}={'dem','rep'};
                mdp(nv,1).label = label;
                mdp(nv,1).T =T;
                mdp(nv,1).a =An;
                mdp(nv,1).B =B;
                mdp(nv,1).C =C;
                mdp(nv,1).D =D;
                mdp(nv,1).E =E;
                %mdp(1,1).beta = precision_belief;
                mdp(nv,1).U =U;
                MDP(nv,1) = spm_MDP_VB_X(mdp(nv,1));
                mdp_array = [mdp_array, nv];

            end

        end

        vote_array = zeros(tv,1);
        dem_count = 0;
        rep_count = 0;

        for nv= 1:tv
           %spm_MDP_VB_trial(MDP(nv,1))
           vote =MDP(nv,1).u(1,1);
           vote_array(nv)=vote;
           if vote == 1
               dem_count = dem_count+1;
           elseif vote ==2
               rep_count = rep_count+1;
           end
        end

        if dem_count>rep_count
            win = 1;
        else
            win = 2;
        end
    winhist(T+1,1)=win;
    disp(state)
    disp(rep_count)
    disp(dem_count)
    victory_margin = [victory_margin,(rep_count - dem_count)];
    
    end








