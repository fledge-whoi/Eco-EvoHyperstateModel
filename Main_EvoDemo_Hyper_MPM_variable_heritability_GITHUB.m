

% Hyperstate model - Quantitative Genetics
clc;clear;close all;

% Three dimensions (m) of stages:
% 1) Age
% 2) Breeding value
% 3) Phenotype
m = 3;

% Set dimensions
w = 2;  % age (or stage)
g = 40; % phenotype (e.g. wing length)

%%  Parameters
%%%%%%%%%%%%%%%%%%
nvital   = 4;    % Number of vital rates (juvenile survival, maturation, adult survival, fertility)
%% Base line life history parameters
%        Fertility | Adult survival | Juvenile survival | Maturation
Data = {[10.9868   , 0.0365         , 0.0965            , 0.9],...   % Species 1
        [4.9958    , 0.2            , 0.25              , 0.572],... % Species 2 
        [0.9987    , 0.8            , 0.385             , 0.4],...   % Species 3  
        [0.3004    , 0.9296         , 0.505             , 0.3],...   % Species 4  
        [0.2286    , 0.95           , 0.8               , 0.07]};    % Species 5
% To guarantee that population is stable with initial life history traits
for idata = 1:5
    SA = Data{idata}(2);
    SJ = Data{idata}(3);
    Y  = Data{idata}(4);
    Data{idata}(1) = (1-SA)*(1-(1-Y)*SJ)/(Y*SJ);
end
ndata = length(Data);
Beta     = 0.15; % Mean effect of phenotype on vital rate
xtheta   = 4;    % Vital rate under selection: 1=f, 2=SJ, 3=SA and 4=Y
xspecies = 5;    % Species characteristics from 1 to 5 

%% Variances: Phenotypic, genetic and environmnetal variance
Vp = 1;       % Phenotypic variance
h2 = 0.2;   % Heritability
Va = h2*Vp;   % Additive genetic variance
Ve = Vp - Va; % Get the environmental variance (Vp = Va + Ve)

%% Create the phenotypic and breeding value bins for the matrix
xp           = 4; % min and max phenotype values
min_pheno    = -xp;
max_pheno    = xp;
min_breeding = -xp;
max_breeding = xp;
mean_pheno   = 0;

% Get the categories of phenotypes
edges_e    = linspace(min_pheno,max_pheno,g+1);
binwidth_e = edges_e(2) - edges_e(1); % Bin width
midpoint_e = edges_e(1:end-1) +  (binwidth_e/2); % Mid phenotype in each bin

% Get the categories of breeding values
edges_a    = edges_e(logical((edges_e<=max_breeding).*(edges_e>=min_breeding)));
binwidth_a = edges_a(2) - edges_a(1); % Bin width
midpoint_a = edges_a(1:end-1) +  (binwidth_a/2); % Mid phenotype in each bin

b = length(midpoint_a); % this is the number of classes for breeding value
midpoint2_a = 2*midpoint_a(1):binwidth_a:2*midpoint_a(end); % I don't get this

%% B matrix
B = speye(b*w*g,b*w*g); % identity matrix

%% P matrix
P = speye(b*w*g,b*w*g); % identity matrix

%% Initial distributions of the population
nind    = 1000; % number of individuals to simulate

N_ve = exp(-(midpoint_e'-midpoint_a).^2./(2*Ve))./sum(exp(-(midpoint_e'-midpoint_a).^2./(2*Ve))*binwidth_e,1);
nbv  = exp(-midpoint_a.^2/(2*Va))./sum(exp(-midpoint_a.^2/(2*Va))*binwidth_a); % Get the proportion of individuals within each class of breeding value
Nbv  = nind*nbv; % Multiply by the total number of individuals in the population to get the number of individuals in each class of breeding value
nz   = N_ve*(nbv')*binwidth_a;
Nz   = nind*nz;
NNinit   = N_ve.*Nbv;

%% Projection time
nit   = 102; %number of years

%% Outputs of the model
M_EBVS = zeros(ndata,nit,nvital); 
SJ_Z   = zeros(ndata,nit,nvital);
V      = zeros(ndata,nvital);
LAMBDA0 = zeros(ndata,nvital);
LAMBDA  = zeros(ndata,nit,nvital);
V_ADAPT = zeros(ndata,nvital);
T_GEN   = zeros(ndata,nvital);
VA = zeros(ndata,nit,nvital);
VP = zeros(ndata,nit,nvital);
H2t = zeros(ndata,nit,nvital);
%% VITAL RATES
for ivital = 1:xtheta

    %% SPECIES
    for idata = 1:xspecies
        F  = Data{idata}(1);
        SA = Data{idata}(2);
        SJ = Data{idata}(3);
        Y  = Data{idata}(4);
        S0=SJ;

        A = [SJ*(1-Y), F; SJ*Y, SA];

        % f=F/SA; % POST BREEDING life cycle
        % f=F/S0; % PRE BREEDING life cycle
        f = F/SJ; % PRE BREEDING life cycle assuming S0=SJ

        [wmat_mean, dmat_mean, vmat_mean] = eig(A);
        [lambda1_mean,imax_mean] = max(diag(dmat_mean)); % this is the dominant eigenvalue, i.e. population growth rate
        LAMBDA0(idata,ivital) = lambda1_mean;
        % Stable age-stage distribution, w = right eigenvector
        w_vec_mean = wmat_mean(:,imax_mean); % right eigenvector corresponding to lambda_1
        w_vec_mean = w_vec_mean./(sum(w_vec_mean)); % normalize w to sum to 1
        v_vec_mean = vmat_mean(:, imax_mean); % left eigenvector corresponding to lambda_1
        v_vec_mean = v_vec_mean./(sum(v_vec_mean)); % normalize v to sum to 1
        S_lambda = v_vec_mean*w_vec_mean'./(v_vec_mean'*w_vec_mean);
        T_GEN(idata,ivital)   = lambda1_mean*(v_vec_mean'*w_vec_mean)/(v_vec_mean'*[0,F;0,0]*w_vec_mean);

        %% Get Theta (parameters)
        S0_pheno = S0;
        SJ_pheno = SJ;
        Y_pheno  = Y;
        SA_pheno = SA;
        f_pheno  = f;
        % F_pheno  = f_pheno*SA_pheno; %POST BREEDING life cycle
        F_pheno  = f_pheno*SJ_pheno;

        if (ivital==1)
            f_pheno  = zeros(b,g);
            %F_pheno  = zeros(b,g); 

            if (idata<=5)
                %% Poisson
                beta = Beta;
                f_pheno  = exp(log(f*ones(b,1)) + beta*midpoint_e);
                dS = Beta*f;
            else
                beta = Beta;
                f_pheno  = invlogit(logit(f*ones(b,1)) + beta*midpoint_e);
                dS = Beta*f*(1-f);
            end
            S_LAMBDA = SJ*S_lambda(1,2);
            F_pheno  =f_pheno*SJ_pheno; %PRE

        elseif(ivital==2)
            SJ_pheno = zeros(b,g);
            beta = Beta;
            SJ_pheno = invlogit(logit(SJ*ones(b,1)) + beta*midpoint_e);
            S_LAMBDA = (1-Y)*S_lambda(1,1)+Y*S_lambda(2,1)+f*S_lambda(1,2);
            dS = Beta*SJ*(1-SJ);
            F_pheno  =f_pheno*SJ_pheno;

        elseif(ivital==3)
            SA_pheno = zeros(b,g);
            beta = Beta;
            SA_pheno = invlogit(logit(SA*ones(b,1)) + beta*midpoint_e);
            %F_pheno=f_pheno.*SA_pheno; %POST BREEDING life cycle
            %S_LAMBDA = S_lambda(2,2)+f*S_lambda(1,2); %POST BREEDING life cycle
            S_LAMBDA = S_lambda(2,2); % PRE BREEDING life cycle
            dS = Beta*SA*(1-SA);

        else
            Y_pheno  = zeros(b,g);
            beta = Beta;
            Y_pheno  = invlogit(logit(Y*ones(b,1)) + beta*midpoint_e);
            S_LAMBDA = SJ*(S_lambda(2,1)-S_lambda(1,1));
            dS = Beta*Y*(1-Y);
        end

        V_ADAPT(idata,ivital) = S_LAMBDA*dS*Va;

        %% Get the matrices U, B, P for Utilde and R, H, M for Ftilde
        %% U matrix
       % Ujk = [SJ_pheno.*(1-Y_pheno),SJ_pheno.*Y_pheno, SA_pheno ];
        idr = [1,2,2]; idc = [1,1,2];
        % Ui = sparse(idr,idc,Ujk);
        nId = length(idr);
        Idr = zeros(nId*b*g,1); Idc = zeros(nId*b*g,1);
        Uu   = zeros(nId*b*g,1);
        for iu = 1:b*g
            Ik = ((iu-1)*nId+1):iu*nId;
            if (ivital ==1)
                Ujk = [SJ_pheno.*(1-Y_pheno),SJ_pheno.*Y_pheno, SA_pheno ];   
            elseif (ivital==2)
                [~,Ptype] = ind2sub([b,g],iu);
                Ujk = [SJ_pheno(1,Ptype).*(1-Y_pheno),SJ_pheno(1,Ptype).*Y_pheno, SA_pheno ];
            elseif(ivital == 3)
                [~,Ptype] = ind2sub([b,g],iu);
                Ujk = [SJ_pheno.*(1-Y_pheno),SJ_pheno.*Y_pheno, SA_pheno(1,Ptype) ];
            elseif(ivital == 4)
                [~,Ptype] = ind2sub([b,g],iu);
                Ujk = [SJ_pheno.*(1-Y_pheno(1,Ptype)),SJ_pheno.*Y_pheno(1,Ptype), SA_pheno ];
            end
            Ui = sparse(idr,idc,Ujk);
            Idr(Ik) = idr+(iu-1)*w; Idc(Ik) = idc+(iu-1)*w;
            Uu(Ik)  = nonzeros(Ui);
        end
        U = sparse(Idr,Idc,Uu,b*w*g,b*w*g);

        %% R matrix
        idr = 1; idc = 2;
        Rjk = F_pheno(1,1);
        Ri = sparse(idr,idc,Rjk);
        nId = length(idr);
        Idr = zeros(nId*b*g,1); Idc = zeros(nId*b*g,1);
        Rr   = zeros(nId*b*g,1);
        for ir = 1:b*g
            Ik = ((ir-1)*nId+1):ir*nId;
            if (ivital==1)||(ivital==2)
            [~,Ptype] = ind2sub([b,g],ir);
            Rjk = F_pheno(1,Ptype);
            Ri = sparse(idr,idc,Rjk);
            end
            Idr(Ik) = idr+(ir-1)*w;Idc(Ik) = idc+(ir-1)*w;
            Rr(Ik) = nonzeros(Ri);
        end
        R = sparse(Idr,Idc,Rr,b*w*g,b*w*g);

        %% H matrix
        % To construct the H matrix, we need the Ga matrix
        Vle = Va;
        GVLE = @(zz) exp(-zz.^2/Va)/sqrt(pi*Va); % Define a new function that gets probabilities from the normal distribution
        % will intervene later in the loop

        %% M matrix
        % To construct the M matrix, we need the Ge matrix
        % For each breeding value, what is the potential distribution of offspring phenotype? - Here we add noise (Ve) to the breeding value to obtain offspring phenotype
        emat = exp(-(midpoint_e'-midpoint_a).^2./(2*Ve))*binwidth_a./sum(exp(-(midpoint_e'-midpoint_a).^2./(2*Ve)).*binwidth_e,1);

        % The matrix M, which is the transmission of maternal phenotype to offspring phenotype goes through the breeding values only
        nId = g*g;
        Idr = zeros(nId*b,1); Idc = zeros(nId*b,1);
        Mm  = zeros(nId*b,1);
        inew = 0;
        for i = 1:b
            iold = inew;
            % Ik = ((i-1)*nId+1):i*nId;
            % Mi = emat(:,i);
            [idr,~,Mi_n] = find(emat(:,i));
            inew =  length(idr);
            idc = reshape(repmat(1:g,inew,1),g*inew,1);
            idr = reshape(repmat(idr,g,1),g*inew,1);
            Ik = ((i-1)*nId+1):(((i-1)*nId+1)+inew*g-1);
            Idr(Ik) = idr+((i-1)*w*g); Idc(Ik) = idc+((i-1)*w*g);
            Mm(Ik) = reshape(repmat(Mi_n,g,1),g*inew,1);
        end
        Idr = nonzeros(Idr);
        Idc = nonzeros(Idc);
        Mm = nonzeros(Mm);
        M = sparse(Idr,Idc,Mm,b*w*g,b*w*g);

        %% Compute N(t)
        %% Initial distribution of the population
        nninit = reshape(permute(NNinit,[2,1]),1,b,g); 
        ninit  = w_vec_mean.*nninit;
        Ninit = ninit(:);

        N      = zeros(w*b*g, nit);
        N(:,1) = Ninit;
        Va_hat = Va;
        Vp_hat = Vp;
        VA(1)  = Va_hat;
        VP(1)  = Vp_hat;
        h2t    = Va_hat/Vp_hat;
        H2t(idata,1,ivital) = h2t;

        for i = 2:nit
            % Live individual transitions
            % Process #1.1 = Transition of live individuals between states
            NU = U*N(:,i-1);

            % Rearrange for process #1.2
            nu = reshape(NU,[w,b,g]);
            nu_r = permute(nu(:,:,:),[2,1,3]); % bwg
            nu_r = nu_r(:);

            % Process #1.2 = Transition of live individuals between breeding value categories (assumed fixed)
            NB = B*nu_r;

            % Rearrange for process #1.3
            nb   = reshape(NB, [b,w,g]);
            nb_r = permute(nb, [3,1,2]);
            nb_r = nb_r(:);

            % Process #1.3 = Transition of live individuals between phenotypic categories (assumed fixed)
            NP = P*nb_r;

            % Rearrange back to initial configuration (stage within breeding value within phenotype)
            np = reshape(NP, [g,b,w]);
            np_r = permute(np, [3,2,1]);
            np_r = np_r(:);
            NU_final = np_r; % Should be the same as NU in the absence of specified transitions between breeding values and phenotypes.

            % Reproduction
            % Process #2.1 = Offspring production
            NF = R*N(:,i-1);

            % Rearrange for process #2.2
            nf = reshape(NF,[w,b,g]);
            n = permute(nf(1,:,:),[3,2,1]); % g,b,1  This is a matrix of offspring produced from maternal phenotypes (rows) and breeding values (columns)
            nn = sum(n,1)./sum(n*binwidth_a,'all');

            % Process #2.2 = Transmission of breeding  value
            Gle = GVLE(0.5*(midpoint2_a));
            H1  = conv2(nn,Gle*binwidth_a);
            H2  = conv2(H1*binwidth_a,n);
            bb2 = (2*midpoint_a(1):(binwidth_a/2):2*midpoint_a(end));
            H3  = (interp1(bb2,H2',midpoint_a));
            nf(1,:,:) = H3*sum(n,'all')./sum(H3,'all');

            % Rearrange for process #2.3
            NH = reshape(permute(nf,[3,1,2]),w*g*b,1);

            % Process # 2.3 = Attribution of offspring phenotype
            %% Compute the variable heritability
            Ve = (1-h2t)*Vp_hat;
            % To construct the M matrix, we need the Ge matrix
            % For each breeding value, what is the potential distribution of offspring phenotype? - Here we add noise (Ve) to the breeding value to obtain offspring phenotype
            emat = exp(-(midpoint_e'-midpoint_a).^2./(2*Ve))*binwidth_a./sum(exp(-(midpoint_e'-midpoint_a).^2./(2*Ve)).*binwidth_e,1);

            % The matrix M, which is the transmission of maternal phenotype to offspring phenotype goes through the breeding values only
            nId = g*g;
            Idr = zeros(nId*b,1); Idc = zeros(nId*b,1);
            Mm  = zeros(nId*b,1);
            inew = 0;
            for im = 1:b
                iold = inew;
                % Ik = ((i-1)*nId+1):i*nId;
                % Mi = emat(:,i);
                [idr,~,Mi_n] = find(emat(:,im));
                inew =  length(idr);
                idc = reshape(repmat(1:g,inew,1),g*inew,1);
                idr = reshape(repmat(idr,g,1),g*inew,1);
                Ik = ((im-1)*nId+1):(((im-1)*nId+1)+inew*g-1);
                Idr(Ik) = idr+((im-1)*w*g); Idc(Ik) = idc+((im-1)*w*g);
                Mm(Ik) = reshape(repmat(Mi_n,g,1),g*inew,1);
            end
            Idr = nonzeros(Idr);
            Idc = nonzeros(Idc);
            Mm = nonzeros(Mm);
            M = sparse(Idr,Idc,Mm,b*w*g,b*w*g);

            NM = M*NH;

            % Rearrange back to initial configuration (state within breeding value within phenotype)
            nm = permute(reshape(NM,[g,w,b]),[2,3,1]);
            nm_r = reshape(nm,w*b*g,1);
            NF_final = nm_r;

            % Add all processes to get total population size
            N_final = NF_final+NU_final;
            N(:,i) = N_final;

            %% Variances computation
            Nn     = reshape(N_final,[w,b,g]); %pop strutured by stage * breding values; phenotype, time
            % Breeding value distribution
            ebvs   = permute(sum(Nn*binwidth_a,[1,3]),[2,1,3]);
            m_ebvs = sum(ebvs.*midpoint_a',1)./sum(ebvs,1);
            % Phenotype distribution
            pheno  = permute(sum(Nn*binwidth_e,[1,2]),[3,1,2]);
            mz     = sum(pheno.*midpoint_e',1)./sum(pheno,1); % phenotype moyen
            % Variances: Additive geentic | Phenotypic
            Va_hat = sum(ebvs.*(midpoint_a'-m_ebvs).^2,1)./sum(ebvs,1); % estimated variance
            VA(idata,i,ivital) = Va_hat;
            Vp_hat = sum(pheno.*(midpoint_e'-mz).^2,1)./sum(pheno,1);   % estimated variance
            VP(idata,i,ivital) = Vp_hat;
            h2t = min(Va_hat/Vp_hat,0.99);
            H2t(idata,i,ivital) = h2t;
        end
        %%
        NN     = reshape(N,[w,b,g,nit]); %pop strutured by stage * breding values; phenotype, time

        ebvs   = permute(sum(NN*binwidth_a,[1,3]),[2,4,1,3]);
        m_ebvs = sum(ebvs.*midpoint_a',1)./sum(ebvs,1);
        m_ebvs2 = ebvs./sum(ebvs.*binwidth_a,1); % distribution of breeding value
        M_EBVS(idata,:,ivital) = m_ebvs;         % mean breeding value 

        ebvsJ   = permute(sum(NN(1,:,:,:)*binwidth_a,[1,3]),[2,4,1,3]);
        m_ebvsJ = sum(ebvsJ.*midpoint_a',1)./sum(ebvsJ,1);
        m_ebvsJ2 = ebvsJ./sum(ebvsJ.*binwidth_a,1); % distribution of breeding value among juveniles
        M_EBVSJ(idata,:,ivital) = m_ebvsJ;          % mean breeding value of juvenile

        ebvsA   = permute(sum(NN(2,:,:,:)*binwidth_a,[1,3]),[2,4,1,3]);
        m_ebvsA = sum(ebvsA.*midpoint_a',1)./sum(ebvsA,1);
        m_ebvsA2 = ebvsA./sum(ebvsA.*binwidth_a,1); % distribution of breeding value among adults
        M_EBVSA(idata,:,ivital) = m_ebvsA;          % mean breeding value of adult

        V(idata,ivital) = mean(m_ebvs(2:nit-2)-m_ebvs(1:nit-2-1));% Ebvs   = cumsum(ebvs*binwidth_a,1);
        
        pheno = permute(sum(NN*binwidth_e,[1,2]),[3,4,1,2]);

        fz   = sum(f_pheno(1,:)'.*pheno,1)./sum(pheno); % mean fertility
        f_Z(idata,:,ivital) = fz;

        SJz   = sum(SJ_pheno(1,:)'.*pheno,1)./sum(pheno); % mean Juvenile survival
        SJ_Z(idata,:,ivital) = SJz;


        SAz   = sum(SA_pheno(1,:)'.*pheno,1)./sum(pheno); % mean adult survival
        SA_Z(idata,:,ivital) = SAz;


        Yz   = sum(Y_pheno(1,:)'.*pheno,1)./sum(pheno); % mean maturation rate
        Y_Z(idata,:,ivital) = Yz;

        mz    = sum(pheno.*midpoint_e',1)./sum(pheno,1); % mean phenotype 

        Va_hat=sum(ebvs.*(midpoint_a'-m_ebvs).^2,1)./sum(ebvs,1); % Genetic variance
        Vp_hat=sum(pheno.*(midpoint_e'-mz).^2,1)./sum(pheno,1);   % Phenotypic variance

    end
    %%
end


[V_ADAPT(idata,ivital) V(idata,ivital)] %theoretical  versus empirical
%%
close all
species_labels = {'species 1', 'species 2', 'species 3', 'species 4', 'species 5'};

colors_rgb =[0.1, 0.4, 0.8;       % Blue
    0.8,0,0;
    1, 165/255, 0;      % Orange
    0.49, 0.18, 0.55     % Purple
    0.5, 0.8, 0.3];       % Dark Green



%%

figure
subplot(2,3,1)
hold on
plot(1:nit,M_EBVS(idata,1:nit,ivital),'color',colors_rgb(idata,:), 'LineWidth',3)
plot(1:nit,M_EBVS(idata,1:nit,ivital)+Va_hat.^0.5,'color',colors_rgb(idata,:), 'LineWidth',1)
plot(1:nit,M_EBVS(idata,1:nit,ivital)-Va_hat.^0.5,'color',colors_rgb(idata,:), 'LineWidth',1)
xlim([0,100])
xlabel('Time (year)','Interpreter','latex','FontSize',16)
ylabel('Mean breeding value','Interpreter','latex','FontSize',16)
%title(['Vital rate is ', num2str(ivital)],'FontSize',16)
set(gca, 'FontSize', 24);

subplot(2,3,2)
plot(midpoint_a, m_ebvs2(:,1),'k','LineWidth',3)
hold on
plot(midpoint_a, m_ebvs2(:,10),'color',[0.49, 0.18, 0.55 ],'LineWidth',3)
plot(midpoint_a, m_ebvs2(:,30),'color',[0.3, 0.7, 0.2],'LineWidth',3)
% plot(midpoint_a, m_ebvs2(:,25),'b','LineWidth',3)
% plot(midpoint_a, m_ebvs2(:,50),'g','LineWidth',3)
plot(midpoint_a, m_ebvs2(:,100),'color',[1, 0.5, 0],'LineWidth',3)
xlim([-xp,xp])
xlabel('Breeding value','Interpreter','latex','FontSize',16)
ylabel('Distribution','Interpreter','latex','FontSize',16)
set(gca, 'FontSize', 24);
%legend('t=0', 't=25','t=50', 't=100')
legend({'t=0', 't=10','t=30', 't=100'}, 'Orientation','horizontal')

subplot(2,3,3)
yyaxis left
plot(1:nit,Va_hat,'color',colors_rgb(idata,:), 'LineWidth',3)
ylabel('Additive genetic variance','Interpreter','latex','FontSize',16)
%ylim([0,0.5])

yyaxis right
plot(1:nit, Vp_hat, ':','color',[1, 0.5, 0], 'LineWidth', 3);
ylabel('Phenotypic variance','Interpreter','latex','FontSize',16)
%ylim([0 1])

xlim([0,100])
xlabel('Time (year)','Interpreter','latex','FontSize',16)
set(gca, 'FontSize', 24);
legend({'Additive genetic variance', 'Phenotypic variance'}, 'Orientation','horizontal')




subplot(2,3,4)
if xtheta==1
    plot(midpoint_e,f_pheno(1,:),'color',colors_rgb(idata,:), 'LineWidth',3)
elseif xtheta==2
    plot(midpoint_e,SJ_pheno(1,:),'color',colors_rgb(idata,:), 'LineWidth',3)
    elseif xtheta==3
    plot(midpoint_e,SA_pheno(1,:),'color',colors_rgb(idata,:), 'LineWidth',3)
else
    plot(midpoint_e,Y_pheno(1,:),'color',colors_rgb(idata,:), 'LineWidth',3)
end
xlim([-xp,xp])
xlabel('Phenotype','Interpreter','latex','FontSize',16)
ylabel('Vital rate','Interpreter','latex','FontSize',16)
set(gca, 'FontSize', 24);

subplot(2,3,5)
if xtheta==1
    plot(f_Z(idata,:,ivital),'color',colors_rgb(idata,:), 'LineWidth',3)
elseif xtheta==2
    plot(SJ_Z(idata,:,ivital),'color',colors_rgb(idata,:), 'LineWidth',3)
    elseif xtheta==3
    plot(SA_Z(idata,:,ivital),'color',colors_rgb(idata,:), 'LineWidth',3)
else
    plot(Y_Z(idata,:,ivital),'color',colors_rgb(idata,:), 'LineWidth',3)
end
xlabel('Time (year)','Interpreter','latex','FontSize',16)
ylabel('Mean vital rate','Interpreter','latex','FontSize',16)
set(gca, 'FontSize', 24);

subplot(2,3,6)
hold on
plot(1:nit,mz , 'color',colors_rgb(idata,:),'LineWidth',3)
plot(1:nit,mz+Vp_hat.^0.5,'color',colors_rgb(idata,:), 'LineWidth',1)
plot(1:nit,mz-Vp_hat.^0.5,'color',colors_rgb(idata,:), 'LineWidth',1)

xlim([0,100])
xlabel('Time (year)','Interpreter','latex','FontSize',16)
ylabel('Mean phenotype','Interpreter','latex','FontSize',16)
%title(['Vital rate is ', num2str(ivital)],'FontSize',16)
set(gca, 'FontSize', 24);