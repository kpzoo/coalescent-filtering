%% Perform Snyder filtering for multi-variable N(t) on real data - no N

% Assumptions and Modifications
% - removed code that removes duplicate sequences
% - added option of dealing with R8s tree
% - removed Nt calculation as it is not optimal and commented plots
% - modified version to work with Garli or other software - assumes tree
% - removed fnid = 7 functionality and ignores Nt posterior
% - starts to use MBE toolbox
% - added confidence bounds and true Nt with last posterior, removed others
% - removed plots including N and lam as no such values with real data
% - modified from batchMV2 to use real fasta sequence data
% - removed sinusoidal uniform sampling code
% - removed minSpace and maxSpace as cells of all fns to single fn arrays
% - the tn time vector is not uniformly sampled and has repetitions
% - added code to calculate a matrix of identifiers for function calling
% - uses x notation and expects certain functional forms
% - used rejection sampling to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clearvars

% Booleans to control functionality and tree 
tic;
plotindiv = 1;
getTree = 0;
linBool = 0; % always 0 for these sims if using snyderRealData

% Boolean to control use of R8s
useR8s = 1;

% Variables for duplicates and mutation clock rate from Pybus2003
remDuplicates = 0; % don't remove duplicates
mu = 0.79*(10^-3);
    
% Functional forms
fnid = 6;
global fnstr;
fnstr = {'x_1e^-x_2t', 'x_1sin(x_2t + x_3) + x_4', 'x_1',...
    'logistic growth', 'piecewise exponential', 'piecewise exponential', 'log piecewise'};
fnstr4 = 'x_1(1 + e^(-x_2x_3))/(1 + e^(-x_2(x_3 - t)) + x_4';
fnstr5 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_4) + x_1e^(-x_2(x_4 - x_3))I(x_4)';
fnstr6 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_3+x_4) + x_1e^(-x_2x_4)I(x_3+x_4)';
numRVset = [2 4 1 4 4 4 4];

%% Extract n and the coalescent times from real fasta data

% Name of newick file
if useR8s
    % Name of R8s newick file
    % Example 63 seqs
    %name = '(((((AF271882i_2153:74.173608,AF271879i_2113:74.173608):10.565146,(AF271878i_2000:6.032860,AF271877i_1999:6.032860):78.705894):208.305988,(AF271875i_1380:47.654265,AF271876i_1797:47.654265):245.390477):22.587139,(AF271887i_3405:78.861740,(AF271880i_2115:2.363943,AF271881i_2116:2.363943):76.497798):236.770140):23.867323,((AF271883i_2432:66.104331,(((AF271886i_3319:38.649464,AF271884i_2659:38.649464):2.893887,AF271885i_3289:41.543351):12.660297,AF271874i_1359:54.203648):11.900683):216.994395,(AF271850i_2673:95.144702,(((AF271828i_0922:72.507660,AF271868i_3471:72.507660):1.516849,(((AF271837i_1801:57.685450,(AF271845i_2386:54.065529,(AF271865i_3463:48.935190,AF271846i_2429:48.935190):5.130339):3.619921):2.127317,(AF271838i_1803:36.973632,AF271835i_1767:36.973632):22.839135,AF271840i_2130:59.812767):12.595671,(((((AF271847i_2438:47.456138,AF271851i_2681:47.456138):11.812831,AF271858i_3318:59.268969):5.857131,(AF271827i_0883:57.962587,AF271860i_3400:57.962587,(AF271857i_2957:55.742694,AF271839i_1997:55.742694):2.219893):7.163512):2.247419,(((((((AF271829i_0923:48.734718,AF271855i_2862:48.734718):2.043564,AF271842i_2141:50.778282):7.315536,(AF271866i_3465:52.002833,(AF271873i_4055:44.128450,AF271843i_2147:44.128450):7.874383):6.090985):1.805479,(AF271849i_2663:53.666824,((AF271826i_0873:42.913002,AF271834i_1339:42.913002):6.647288,AF271856i_2926:49.560290):4.106534,(AF271833i_1240:44.619278,(AF271870i_4033:14.776011,AF271844i_2150:14.776011,AF271861i_3452:14.776011):29.843267):9.047547):6.232472,AF271825i_0800:59.899297):2.301461,((AF271852i_2685:47.760340,AF271853i_2852:47.760340):4.114632,(AF271867i_3468:48.729916,(AF271862i_3458:26.684969,AF271864i_3461:26.684969):22.044947,AF271859i_3393:48.729916):3.145057):10.325785):2.519314,((AF271832i_1226:7.279427,AF271836i_1796:7.279427):56.766598,AF271872i_4053:64.046024):0.674047):1.208066,(((AF271869i_4020:13.782170,AF271871i_4036:13.782170):31.500325,AF271830i_1150:45.282495):8.541679,AF271863i_3460:53.824175):12.103962,(AF271831i_1164:58.249836,AF271841i_2134:58.249836):7.678301):1.445382):4.100041,AF271848i_2446:71.473560):0.934877):1.616071):6.684342,AF271854i_2856:80.708851):14.435852):187.954023):56.400479)root;';
    
    % Best 63 sequence of 50 garli trees
    % This one isn't so good, constrained 100-500
    %name = '(((((AF271882i_2153:21.868834,AF271879i_2113:21.868834):3.309316,(AF271878i_2000:1.726499,AF271877i_1999:1.726499):23.451652):61.563240,(AF271875i_1380:14.044640,AF271876i_1797:14.044640):72.696751):6.478811,(AF271887i_3405:23.408771,(AF271880i_2115:0.651426,AF271881i_2116:0.651426):22.757345):69.811430):6.779909,((AF271883i_2432:19.361451,(((AF271886i_3319:11.370143,AF271884i_2659:11.370143):1.325783,AF271885i_3289:12.695926):3.320500,AF271874i_1359:16.016426):3.345025):63.692615,(AF271850i_2673:27.891118,((AF271828i_0922:21.690980,AF271868i_3471:21.690980,(((AF271837i_1801:17.086822,(AF271845i_2386:16.028913,(AF271865i_3463:14.475485,AF271846i_2429:14.475485):1.553428):1.057909):0.540257,(AF271838i_1803:10.917094,AF271835i_1767:10.917094):6.709985,AF271840i_2130:17.627079):3.517578,(((((AF271847i_2438:13.733489,AF271851i_2681:13.733489):3.466612,AF271858i_3318:17.200101):1.763803,(AF271827i_0883:16.771621,AF271860i_3400:16.771621,(AF271857i_2957:15.882001,AF271839i_1997:15.882001):0.889620):2.192283):0.806938,(((((((AF271829i_0923:14.430932,AF271855i_2862:14.430932):0.815614,AF271842i_2141:15.246546):1.986258,(AF271866i_3465:15.464985,(AF271873i_4055:13.118136,AF271843i_2147:13.118136):2.346849):1.767819):0.389607,(AF271849i_2663:15.963346,((AF271826i_0873:12.650436,AF271834i_1339:12.650436):2.002493,AF271856i_2926:14.652929):1.310417,(AF271833i_1240:13.303260,(AF271870i_4033:4.383157,AF271844i_2150:4.383157,AF271861i_3452:4.383157):8.920103):2.660086):1.659064,AF271825i_0800:17.622410):0.660183,((AF271852i_2685:14.292180,AF271853i_2852:14.292180):0.961183,(AF271867i_3468:14.432707,(AF271862i_3458:7.721612,AF271864i_3461:7.721612):6.711095,AF271859i_3393:14.432707):0.820656):3.029230):0.726944,(AF271832i_1226:2.043016,AF271836i_1796:2.043016):16.966521,AF271872i_4053:19.009537):0.355463,(((AF271869i_4020:4.068606,AF271871i_4036:4.068606):9.249341,AF271830i_1150:13.317947):2.230077,AF271863i_3460:15.548025):3.816975,(AF271831i_1164:17.225120,AF271841i_2134:17.225120):2.139880):0.405843):1.088750,AF271848i_2446:20.859593):0.285064):0.546323):1.960792,AF271854i_2856:23.651772):4.239346):55.162948):16.946045)root;';
    % This one is better - constrained 248-322 from Pybus
    %name = '(((((AF271882i_2153:68.026652,AF271879i_2113:68.026652):10.294180,(AF271878i_2000:5.370564,AF271877i_1999:5.370564):72.950269):191.502717,(AF271875i_1380:43.688195,AF271876i_1797:43.688195):226.135354):20.153419,(AF271887i_3405:72.816884,(AF271880i_2115:2.026370,AF271881i_2116:2.026370):70.790515):217.160084):21.090035,((AF271883i_2432:60.227016,(((AF271886i_3319:35.368723,AF271884i_2659:35.368723):4.124068,AF271885i_3289:39.492792):10.328967,AF271874i_1359:49.821759):10.405258):198.126489,(AF271850i_2673:86.759967,((AF271828i_0922:67.473406,AF271868i_3471:67.473406,(((AF271837i_1801:53.151404,(AF271845i_2386:49.860603,(AF271865i_3463:45.028407,AF271846i_2429:45.028407):4.832196):3.290801):1.680561,(AF271838i_1803:33.959438,AF271835i_1767:33.959438):20.872527,AF271840i_2130:54.831964):10.942013,(((((AF271847i_2438:42.720305,AF271851i_2681:42.720305):10.783475,AF271858i_3318:53.503780):5.486603,(AF271827i_0883:52.170921,AF271860i_3400:52.170921,(AF271857i_2957:49.403609,AF271839i_1997:49.403609):2.767312):6.819462):2.510116,(((((((AF271829i_0923:44.889818,AF271855i_2862:44.889818):2.537102,AF271842i_2141:47.426920):6.178587,(AF271866i_3465:48.106412,(AF271873i_4055:40.806147,AF271843i_2147:40.806147):7.300265):5.499095):1.211936,(AF271849i_2663:49.656647,((AF271826i_0873:39.351288,AF271834i_1339:39.351288):6.229088,AF271856i_2926:45.580375):4.076272,(AF271833i_1240:41.382007,(AF271870i_4033:13.634539,AF271844i_2150:13.634539,AF271861i_3452:13.634539):27.747468):8.274640):5.160796,AF271825i_0800:54.817442):2.053609,((AF271852i_2685:44.458207,AF271853i_2852:44.458207):2.989920,(AF271867i_3468:44.895339,(AF271862i_3458:24.019359,AF271864i_3461:24.019359):20.875979,AF271859i_3393:44.895339):2.552789):9.422924):2.261280,(AF271832i_1226:6.355141,AF271836i_1796:6.355141):52.777190,AF271872i_4053:59.132331):1.105726,(((AF271869i_4020:12.656076,AF271871i_4036:12.656076):28.771616,AF271830i_1150:41.427693):6.937028,AF271863i_3460:48.364720):11.873336,(AF271831i_1164:53.581604,AF271841i_2134:53.581604):6.656452):1.262443):3.386738,AF271848i_2446:64.887237):0.886741):1.699428):6.099370,AF271854i_2856:73.572776):13.187191):171.593538):52.713498)root;';
    
    % Constrained 248-322 but nrates = 3
    name = '(((((AF271882i_2153:62.364393,AF271879i_2113:62.364393):9.436276,(AF271878i_2000:4.924179,AF271877i_1999:4.924179):66.876490):175.602359,(AF271875i_1380:40.058496,AF271876i_1797:40.058496):207.344531):18.479688,(AF271887i_3405:66.767417,(AF271880i_2115:1.858035,AF271881i_2116:1.858035):64.909382):199.115298):19.338876,((AF271883i_2432:55.226772,(((AF271886i_3319:32.432086,AF271884i_2659:32.432086):3.781891,AF271885i_3289:36.213978):9.471452,AF271874i_1359:45.685430):9.541342):181.662832,(AF271850i_2673:79.552491,((AF271828i_0922:61.868688,AF271868i_3471:61.868688,(((AF271837i_1801:48.741694,(AF271845i_2386:45.723781,(AF271865i_3463:41.292123,AF271846i_2429:41.292123):4.431658):3.017913):1.540965,(AF271838i_1803:31.141236,AF271835i_1767:31.141236):19.141423,AF271840i_2130:50.282659):10.028297,(((((AF271847i_2438:39.173051,AF271851i_2681:39.173051):9.888203,AF271858i_3318:49.061253):5.030900,(AF271827i_0883:47.837898,AF271860i_3400:47.837898,(AF271857i_2957:45.300396,AF271839i_1997:45.300396):2.537502):6.254255):2.301486,(((((((AF271829i_0923:41.161586,AF271855i_2862:41.161586):2.326417,AF271842i_2141:43.488003):5.665451,(AF271866i_3465:44.110995,(AF271873i_4055:37.417001,AF271843i_2147:37.417001):6.693994):5.042459):1.111321,(AF271849i_2663:45.532553,((AF271826i_0873:36.083004,AF271834i_1339:36.083004):5.711790,AF271856i_2926:41.794794):3.737759,(AF271833i_1240:37.945029,(AF271870i_4033:12.502058,AF271844i_2150:12.502058,AF271861i_3452:12.502058):25.442970):7.587524):4.732222,AF271825i_0800:50.264775):1.883254,((AF271852i_2685:40.766276,AF271853i_2852:40.766276):2.741707,(AF271867i_3468:41.167061,(AF271862i_3458:22.024539,AF271864i_3461:22.024539):19.142522,AF271859i_3393:41.167061):2.340922):8.640046):2.073820,(AF271832i_1226:5.827327,AF271836i_1796:5.827327):48.394522,AF271872i_4053:54.221849):1.014109,(((AF271869i_4020:11.604912,AF271871i_4036:11.604912):26.382202,AF271830i_1150:37.987113):6.361073,AF271863i_3460:44.348187):10.887772,(AF271831i_1164:49.132127,AF271841i_2134:49.132127):6.103831):1.157681):3.104557,AF271848i_2446:59.498196):0.812760):1.557732):5.592167,AF271854i_2856:67.460856):12.091635):157.337113):48.331987)root;';
    
    % Best 68 sequence of 50 garli trees, constrained 300-800
    %name = '((AF271876i_1797:277.582528,(AF271875i_1380:255.414011,(((AF271878i_2000:4.090012,AF271877i_1999:4.090012):50.868515,(AF271882i_2153:48.713581,AF271879i_2113:48.713581):6.244946):115.191608,(((AF271883i_2432:41.740881,(((AF271884i_2659:25.664530,AF271886i_3319:25.664530):2.941489,AF271885i_3289:28.606019):6.637878,AF271874i_1359:35.243897):6.496984):93.411928,(AF271850i_2673:60.690948,(AF271854i_2856:53.013209,(AF271845i_2386:51.289737,(AF271828i_0922:50.121516,((((AF271840i_2130:36.851760,(AF271835i_1767:23.169957,AF271838i_1803:23.169957):13.681804,AF271837i_1801:36.851760):3.366663,(AF271865i_3463:34.513158,AF271846i_2429:34.513158):5.705265):4.978098,(AF271848i_2446:42.914172,AF271868i_3471:42.914172):2.282350):3.820189,((((AF271851i_2681:33.900047,AF271847i_2438:33.900047):8.993098,AF271858i_3318:42.893145):3.964366,((AF271832i_1226:5.061813,AF271836i_1796:5.061813):40.379303,((((((AF271852i_2685:32.010757,(AF271862i_3458:17.815318,AF271864i_3461:17.815318):14.195440):2.420026,(AF271859i_3393:32.316431,AF271853i_2852:32.316431):2.114353):2.710568,AF271867i_3468:37.141351):4.246899,(AF271825i_0800:40.078788,((AF271866i_3465:35.271022,(AF271843i_2147:29.948887,AF271873i_4055:29.948887):5.322135):3.984078,(AF271842i_2141:35.056108,(AF271829i_0923:33.194410,AF271855i_2862:33.194410):1.861698):4.198991):0.823689,(AF271849i_2663:36.544298,(AF271856i_2926:33.577160,(AF271826i_0873:28.998557,AF271834i_1339:28.998557):4.578602):2.967138,(AF271833i_1240:30.517346,(AF271870i_4033:10.656787,AF271861i_3452:10.656787,AF271844i_2150:10.656787):19.860559):6.026952):3.534490):1.309462):2.025622,(AF271872i_4053:38.057396,AF271831i_1164:38.057396):5.356475,(AF271857i_2957:35.681669,(AF271839i_1997:33.030245,AF271860i_3400:33.030245):2.651425):7.732202):0.656205,AF271841i_2134:44.070077,(((AF271869i_4020:10.226527,AF271871i_4036:10.226527):21.710871,AF271830i_1150:31.937398):5.707849,AF271863i_3460:37.645247):6.424829):1.371039):1.416395):1.050545,AF271827i_0883:47.908056):1.108654):1.104806):1.168221):1.723472):7.677739):74.461860):26.542651,(AF271887i_3405:46.459032,(AF271881i_2116:1.445447,AF271880i_2115:1.445447):45.013586):115.236428):8.454675):85.263876):22.168516):406.452263,(AF271820i_1382:75.435446,((AF271822i_2152:39.555665,(AF271821i_2004:31.770574,AF271824i_3664:31.770574):7.785090):3.479310,AF271823i_2440:43.034975):32.400471):608.599345)root;';
    % This one is constrained 512-735 from Pybus
    %name = '((AF271876i_1797:298.026214,(AF271875i_1380:274.225009,(((AF271878i_2000:4.391238,AF271877i_1999:4.391238):54.614935,(AF271882i_2153:52.301292,AF271879i_2113:52.301292):6.704881):123.675360,(((AF271883i_2432:44.815055,(((AF271884i_2659:27.554697,AF271886i_3319:27.554697):3.158127,AF271885i_3289:30.712824):7.126751,AF271874i_1359:37.839574):6.975481):100.291630,(AF271850i_2673:65.160776,(AF271854i_2856:56.917578,(AF271845i_2386:55.067174,(AF271828i_0922:53.812915,((((AF271840i_2130:39.565855,(AF271835i_1767:24.876401,AF271838i_1803:24.876401):14.689454,AF271837i_1801:39.565855):3.614614,(AF271865i_3463:37.055017,AF271846i_2429:37.055017):6.125452):5.344730,(AF271848i_2446:46.074757,AF271868i_3471:46.074757):2.450442):4.101542,((((AF271851i_2681:36.396751,AF271847i_2438:36.396751):9.655431,AF271858i_3318:46.052182):4.256338,((AF271832i_1226:5.434610,AF271836i_1796:5.434610):43.353198,((((((AF271852i_2685:34.368318,(AF271862i_3458:19.127399,AF271864i_3461:19.127399):15.240919):2.598258,(AF271859i_3393:34.696503,AF271853i_2852:34.696503):2.270073):2.910198,AF271867i_3468:39.876774):4.559678,(AF271825i_0800:43.030550,((AF271866i_3465:37.868697,(AF271843i_2147:32.154592,AF271873i_4055:32.154592):5.714105):4.277501,(AF271842i_2141:37.637955,(AF271829i_0923:35.639145,AF271855i_2862:35.639145):1.998810):4.508242):0.884353,(AF271849i_2663:39.235749,(AF271856i_2926:36.050084,(AF271826i_0873:31.134272,AF271834i_1339:31.134272):4.915812):3.185665,(AF271833i_1240:32.764918,(AF271870i_4033:11.441649,AF271861i_3452:11.441649,AF271844i_2150:11.441649):21.323269):6.470831):3.794801):1.405902):2.174806,(AF271872i_4053:40.860285,AF271831i_1164:40.860285):5.750974,(AF271857i_2957:38.309588,(AF271839i_1997:35.462889,AF271860i_3400:35.462889):2.846699):8.301671):0.704534,AF271841i_2134:47.315793,(((AF271869i_4020:10.979700,AF271871i_4036:10.979700):23.309855,AF271830i_1150:34.289555):6.128227,AF271863i_3460:40.417781):6.898012):1.472015):1.520711):1.127917,AF271827i_0883:51.436436):1.190305):1.186174):1.254259):1.850404):8.243198):79.945909):28.497494,(AF271887i_3405:49.880698,(AF271881i_2116:1.551902,AF271880i_2115:1.551902):48.328796):123.723481):9.077354):91.543476):23.801205):436.387081,(AF271820i_1382:80.991194,((AF271822i_2152:42.468902,(AF271821i_2004:34.110447,AF271824i_3664:34.110447):8.358455):3.735557,AF271823i_2440:46.204459):34.786735):653.422101)root;';

    % Constrained 512-735 but nrates = 3
    %name = '((AF271876i_1797:279.848520,(AF271875i_1380:257.497223,(((AF271878i_2000:4.122927,AF271877i_1999:4.122927):51.280671,(AF271882i_2153:49.107873,AF271879i_2113:49.107873):6.295725):116.131882,(((AF271883i_2432:42.075914,(((AF271884i_2659:25.870245,AF271886i_3319:25.870245):2.965025,AF271885i_3289:28.835270):6.691293,AF271874i_1359:35.526563):6.549350):94.170468,(AF271850i_2673:61.167806,(AF271854i_2856:53.427596,(AF271845i_2386:51.690526,(AF271828i_0922:50.513333,((((AF271840i_2130:37.143247,(AF271835i_1767:23.353673,AF271838i_1803:23.353673):13.789574,AF271837i_1801:37.143247):3.392520,(AF271865i_3463:34.786052,AF271846i_2429:34.786052):5.749715):5.015967,(AF271848i_2446:43.251846,AF271868i_3471:43.251846):2.299888):3.848643,((((AF271851i_2681:34.167169,AF271847i_2438:34.167169):9.063410,AF271858i_3318:43.230579):3.995091,((AF271832i_1226:5.102081,AF271836i_1796:5.102081):40.697438,((((((AF271852i_2685:32.264032,(AF271862i_3458:17.956914,AF271864i_3461:17.956914):14.307119):2.438777,(AF271859i_3393:32.572004,AF271853i_2852:32.572004):2.130805):2.732922,AF271867i_3468:37.435732):4.282244,(AF271825i_0800:40.398817,((AF271866i_3465:35.552933,(AF271843i_2147:30.188350,AF271873i_4055:30.188350):5.364583):4.015747,(AF271842i_2141:35.336422,(AF271829i_0923:33.459884,AF271855i_2862:33.459884):1.876538):4.232258):0.830138,(AF271849i_2663:36.836493,(AF271856i_2926:33.845754,(AF271826i_0873:29.230598,AF271834i_1339:29.230598):4.615156):2.990739,(AF271833i_1240:30.761501,(AF271870i_4033:10.742179,AF271861i_3452:10.742179,AF271844i_2150:10.742179):20.019322):6.074991):3.562325):1.319159):2.040221,(AF271872i_4053:38.359604,AF271831i_1164:38.359604):5.398594,(AF271857i_2957:35.965229,(AF271839i_1997:33.292860,AF271860i_3400:33.292860):2.672370):7.792968):0.660782,AF271841i_2134:44.418979,(((AF271869i_4020:10.308232,AF271871i_4036:10.308232):21.883468,AF271830i_1150:32.191700):5.752708,AF271863i_3460:37.944408):6.474572):1.380539):1.426152):1.058024,AF271827i_0883:48.283694):1.116683):1.112957):1.177192):1.737071):7.740210):75.078575):26.763893,(AF271887i_3405:46.834915,(AF271881i_2116:1.457106,AF271880i_2115:1.457106):45.377808):116.175360):8.525206):85.961743):22.351297):409.693370,(AF271820i_1382:76.045210,((AF271822i_2152:39.876436,(AF271821i_2004:32.027942,AF271824i_3664:32.027942):7.848495):3.507660,AF271823i_2440:43.384096):32.661114):613.496681)root;';
    % Fixed at mean of 619
    %name = '((AF271876i_1797:251.191295,(AF271875i_1380:231.130457,(((AF271878i_2000:3.701153,AF271877i_1999:3.701153):46.032178,(AF271882i_2153:44.082125,AF271879i_2113:44.082125):5.651206):104.239735,(((AF271883i_2432:37.772355,(((AF271884i_2659:23.224468,AF271886i_3319:23.224468):2.661826,AF271885i_3289:25.886294):6.006780,AF271874i_1359:31.893074):5.879281):84.530765,(AF271850i_2673:54.920740,(AF271854i_2856:47.972963,(AF271845i_2386:46.413351,(AF271828i_0922:45.356199,((((AF271840i_2130:33.348069,(AF271835i_1767:20.967067,AF271838i_1803:20.967067):12.381002,AF271837i_1801:33.348069):3.046577,(AF271865i_3463:31.231810,AF271846i_2429:31.231810):5.162835):4.504804,(AF271848i_2446:38.834095,AF271868i_3471:38.834095):2.065355):3.456983,((((AF271851i_2681:30.676990,AF271847i_2438:30.676990):8.138077,AF271858i_3318:38.815067):3.587453,((AF271832i_1226:4.580559,AF271836i_1796:4.580559):36.540230,((((((AF271852i_2685:28.967326,(AF271862i_3458:16.121522,AF271864i_3461:16.121522):12.845804):2.189941,(AF271859i_3393:29.243937,AF271853i_2852:29.243937):1.913330):2.452860,AF271867i_3468:33.610127):3.843123,(AF271825i_0800:36.268286,((AF271866i_3465:31.917619,(AF271843i_2147:27.101488,AF271873i_4055:27.101488):4.816131):3.605290,(AF271842i_2141:31.723139,(AF271829i_0923:30.038442,AF271855i_2862:30.038442):1.684697):3.799771):0.745377,(AF271849i_2663:33.069839,(AF271856i_2926:30.384802,(AF271826i_0873:26.241511,AF271834i_1339:26.241511):4.143291):2.685037,(AF271833i_1240:27.615901,(AF271870i_4033:9.643590,AF271861i_3452:9.643590,AF271844i_2150:9.643590):17.972311):5.453938):3.198447):1.184964):1.833035,(AF271872i_4053:34.439079,AF271831i_1164:34.439079):4.847207,(AF271857i_2957:32.289226,(AF271839i_1997:29.889886,AF271860i_3400:29.889886):2.399339):6.997060):0.593816,AF271841i_2134:39.880102,(((AF271869i_4020:9.254237,AF271871i_4036:9.254237):19.646704,AF271830i_1150:28.900941):5.165174,AF271863i_3460:34.066115):5.813987):1.240687):1.281731):0.950664,AF271827i_0883:43.353184):1.003248):0.999766):1.057152):1.559612):6.947776):67.382380):24.019102,(AF271887i_3405:42.041927,(AF271881i_2116:1.308020,AF271880i_2115:1.308020):40.733907):104.280295):7.650844):77.157390):20.060839):367.808705,(AF271820i_1382:68.263401,((AF271822i_2152:35.794900,(AF271821i_2004:28.749979,AF271824i_3664:28.749979):7.044921):3.148513,AF271823i_2440:38.943413):29.319988):550.736599)root;';
    
    % Single rate LF models - 63 constrained to 248-322
    %name = '(((((AF271882i_2153:54.264659,AF271879i_2113:54.264659):8.211637,(AF271878i_2000:4.284082,AF271877i_1999:4.284082):58.192214):152.761154,(AF271875i_1380:34.849942,AF271876i_1797:34.849942):180.387508):16.076323,(AF271887i_3405:58.085814,(AF271880i_2115:1.616429,AF271881i_2116:1.616429):56.469385):173.227959):16.823459,((AF271883i_2432:48.042915,(((AF271886i_3319:28.213527,AF271884i_2659:28.213527):3.289757,AF271885i_3289:31.503284):8.239387,AF271874i_1359:39.742671):8.300244):158.044915,(AF271850i_2673:69.208169,((AF271828i_0922:53.823336,AF271868i_3471:53.823336,(((AF271837i_1801:42.398718,(AF271845i_2386:39.773656,(AF271865i_3463:35.919027,AF271846i_2429:35.919027):3.854629):2.625062):1.340578,(AF271838i_1803:27.089343,AF271835i_1767:27.089343):16.649953,AF271840i_2130:43.739296):8.728412,(((((AF271847i_2438:34.077861,AF271851i_2681:34.077861):8.601946,AF271858i_3318:42.679807):4.376647,(AF271827i_0883:41.616591,AF271860i_3400:41.616591,(AF271857i_2957:39.409114,AF271839i_1997:39.409114):2.207476):5.439864):2.002312,(((((((AF271829i_0923:35.808475,AF271855i_2862:35.808475):2.023839,AF271842i_2141:37.832313):4.928640,(AF271866i_3465:38.374342,(AF271873i_4055:32.550942,AF271843i_2147:32.550942):5.823400):4.386612):0.966758,(AF271849i_2663:39.610960,((AF271826i_0873:31.390406,AF271834i_1339:31.390406):4.968925,AF271856i_2926:36.359330):3.251630,(AF271833i_1240:33.010304,(AF271870i_4033:10.876231,AF271844i_2150:10.876231,AF271861i_3452:10.876231):22.134073):6.600656):4.116751,AF271825i_0800:43.727711):1.638157,((AF271852i_2685:35.464180,AF271853i_2852:35.464180):2.385050,(AF271867i_3468:35.812878,(AF271862i_3458:19.160171,AF271864i_3461:19.160171):16.652707,AF271859i_3393:35.812878):2.036352):7.516639):1.803816,(AF271832i_1226:5.069477,AF271836i_1796:5.069477):42.100208,AF271872i_4053:47.169685):0.882034,(((AF271869i_4020:10.095715,AF271871i_4036:10.095715):22.951033,AF271830i_1150:33.046748):5.533646,AF271863i_3460:38.580394):9.471325,(AF271831i_1164:42.741887,AF271841i_2134:42.741887):5.309832):1.007047):2.701591,AF271848i_2446:51.760357):0.707350):1.355628):4.865449,AF271854i_2856:58.688786):10.519383):136.879661):42.049403)root;';
    
    % Single rate LF models - 68 constrained to 512-735
    %name = '((AF271876i_1797:207.771543,(AF271875i_1380:191.178327,(((AF271878i_2000:3.061389,AF271877i_1999:3.061389):38.075271,(AF271882i_2153:36.462295,AF271879i_2113:36.462295):4.674365):86.221342,(((AF271883i_2432:31.243202,(((AF271884i_2659:19.209995,AF271886i_3319:19.209995):2.201715,AF271885i_3289:21.411710):4.968476,AF271874i_1359:26.380186):4.863016):69.919172,(AF271850i_2673:45.427398,(AF271854i_2856:39.680582,(AF271845i_2386:38.390557,(AF271828i_0922:37.516139,((((AF271840i_2130:27.583678,(AF271835i_1767:17.342798,AF271838i_1803:17.342798):10.240880,AF271837i_1801:27.583678):2.519960,(AF271865i_3463:25.833226,AF271846i_2429:25.833226):4.270412):3.726125,(AF271848i_2446:32.121416,AF271868i_3471:32.121416):1.708347):2.859425,((((AF271851i_2681:25.374310,AF271847i_2438:25.374310):6.731367,AF271858i_3318:32.105677):2.967343,((AF271832i_1226:3.788785,AF271836i_1796:3.788785):30.224057,((((((AF271852i_2685:23.960169,(AF271862i_3458:13.334831,AF271864i_3461:13.334831):10.625339):1.811398,(AF271859i_3393:24.188967,AF271853i_2852:24.188967):1.582601):2.028870,AF271867i_3468:27.800438):3.178819,(AF271825i_0800:29.999120,((AF271866i_3465:26.400489,(AF271843i_2147:22.416852,AF271873i_4055:22.416852):3.983638):2.982097,(AF271842i_2141:26.239625,(AF271829i_0923:24.846138,AF271855i_2862:24.846138):1.393488):3.142960):0.616534,(AF271849i_2663:27.353541,(AF271856i_2926:25.132627,(AF271826i_0873:21.705526,AF271834i_1339:21.705526):3.427101):2.220914,(AF271833i_1240:22.842346,(AF271870i_4033:7.976644,AF271861i_3452:7.976644,AF271844i_2150:7.976644):14.865701):4.511196):2.645579):0.980137):1.516185,(AF271872i_4053:28.486101,AF271831i_1164:28.486101):4.009341,(AF271857i_2957:26.707861,(AF271839i_1997:24.723260,AF271860i_3400:24.723260):1.984601):5.787581):0.491172,AF271841i_2134:32.986615,(((AF271869i_4020:7.654593,AF271871i_4036:7.654593):16.250667,AF271830i_1150:23.905260):4.272346,AF271863i_3460:28.177606):4.809008):1.026228):1.060177):0.786337,AF271827i_0883:35.859357):0.829832):0.826951):0.874418):1.290025):5.746816):55.734976):19.867272,(AF271887i_3405:34.774757,(AF271881i_2116:1.081922,AF271880i_2115:1.081922):33.692835):86.254890):6.328355):63.820324):16.593216):304.231014,(AF271820i_1382:56.463709,((AF271822i_2152:29.607561,(AF271821i_2004:23.780392,AF271824i_3664:23.780392):5.827169):2.604276,AF271823i_2440:32.211838):24.251872):455.538847)root;';
else
    % Name of Garli file
    name = '63best50';
end

% Read tree, assume tree in evolutionary distances of some kind
treeSeq = phytreeread(name);

% Get the branch to leaf distances and the no. lineages
[M, ~, D] = getmatrix(treeSeq);
n = get(treeSeq, 'NumLeaves');
nold = n;

% Ensure only 63 or 68 sequences
if ~any(nold == [63 68])
    error('Incorrect data supplied');
end

% Get true branch lengths using the connectivity matrix and distances
[T, dcomp] = getCoalDists(M, n, D);

% Obtain full coalescent time vector - converting distance to times
evDist = [0 T];
evDist = sort(evDist);
if useR8s
    % Branches already in time
    tcoal = evDist;
else
    % Branches in number of substitutions per site per year
    tcoal = evDist/mu;
end
twait = diff(tcoal);

% Other way to account for coalescent duplicates instead of maintaining
% the filter coalescent
if remDuplicates
    if(any(diff(tcoal) == 0))
        tcoal = unique(tcoal);
        n = length(tcoal);
        disp(['Duplicate coalescents means n goes from ' num2str(nold) ' to ' num2str(n)]);
    end
end

%% Parameter settings

% Set number of parameters, dimensions (mi) and data length
mi = [20 20 20 20];
numRV = length(mi);
nData = n-1;
m = prod(mi);

% Check for pre-specified parameters not to be estimated - mi values of 1
specParam = find(mi == 1);
nEstRV = numRV - length(specParam);

% Get coalescent binomial factors for lineages and change notation
nset = n:-1:2;
fac = nset.*(nset-1)/2;

% Space sets for functions
switch(fnid)
    case 1
        % Exponential
        minSpace = [100*n 0.1];
        maxSpace = [1000*n 10];
    case 2
        % Sinusoidal
        minSpace = [100*n 0.1 0 1100*n];
        maxSpace = [1000*n 10 pi/2 1200*n];
    case 3
        % Constant
        minSpace = 100*n;
        maxSpace = 1000*n;
    case 4
        % Logistic
        minSpace = [1000 0 50 50];
        maxSpace = [10000 10 max(tcoal) 500];
    case 5
        % Pybus piecewise exponential
        minSpace = [5000 0 0 0];
        maxSpace = [15000 1 200 200];
    case 6
        % Pybus piecewise exponential but specify x and y-x vs y
        minSpace = [1000 0 0 0];
        maxSpace = [20000 0.75 80 100];
end

% Modify spaces if specify parameters so that they are set to min values
if nEstRV > 0
    minSpace(specParam) = maxSpace(specParam);
end

% Function identifier structure
fn.id = fnid;
fn.name = fnstr(fnid);
fn.nData = nData;
fn.mi = mi;
fn.m = m;
fn.fac = fac;
fn.numRV = numRV;
fn.nEstRV = nEstRV;

%% Initialise and define key parameters and storage variables

% Check numRV compatible with fntype
if numRV ~= numRVset(fnid)
    error('Initial settings on RV space inconsistent');
end

% Set uniform prior q0
q0 = ones(1, m)/m;

% Get the random variable sets
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
end
fn.xset = xset;

% Create a matrix of identifiers to tell which xset{i} values are used for
% each entry in N(t) and lam(t) calculations
IDMx = zeros(numRV, m);
% Initialise with first variable which has no element repetitions
idxset = 1:mi(1);
IDMx(1, :) = repmat(idxset, 1, m/mi(1));
for i = 2:numRV
    % For further variables numReps gives the number of set repetitions
    % while kronVec gives the number of element repetitions
    idxset = 1:mi(i);
    numReps = m/prod(mi(1:i));
    kronVec = ones(1, prod(mi(1:i-1)));
    IDMx(i, :) = repmat(kron(idxset, kronVec), 1, numReps);
end

% Get the values corresponding to the matrix
xsetMx = zeros(numRV, m);
for i = 1:numRV
    xsetMx(i, :) = xset{i}(IDMx(i, :));
end
fn.IDMx = IDMx;
fn.xsetMx = xsetMx;

%% Filter the coalescent times from the fasta data

% Main function to perform Snyder estimation and return MMSE estimates
[xhat, Nhat, lamhat, tn, nLin, statsStruc, qn] = snyderRealDataDuplicate(fn, q0, tcoal);
tnmin = min(tn);
tnmax = max(tn);
tnlen = length(tn);
paramhat = xhat(end, :)';
fn.paramhat = paramhat;

% Extract outputs
xhatmean = statsStruc.xhatmean;
xstdlb = statsStruc.xstdlb;
xstdub = statsStruc.xstdub;

% Marginalise the posterior across sim times - take 20 of them
nPoster = 20;
idpost = round(linspace(1, length(tn), nPoster));
qmarg = cell(1, nPoster);
probSums = cell(1, nPoster);
for i = 1:nPoster
    [qmarg{i}, probSums{i}] = marginalise(numRV, IDMx, qn(idpost(i), :), mi);
end
qnlast = qn(end, :);
clear qn

% Reorganise into probability matrices for each variable
qnp = cell(1, numRV);
for i = 1:numRV
    qnp{i} = zeros(nPoster, mi(i));
    for j = 1:nPoster
        qnp{i}(j, :) = qmarg{j}{i};
    end
end

% Get optimal values backward in time like coalescent - it is a smooth
% function based on the final posterior
t = linspace(tcoal(1), tcoal(end), 1000);
[Nf2, Nbnd2] = getPopBackward(t, qnlast, fn);

% Reproduce Pybus2003 plot with time reorder and confidence intervals
if fnid == 5 || 6
    % Get other parameters of interest
    pm = [paramhat; 0];
    if fnid == 6
        % Param 4 is y-x
        pm(5) = paramhat(1)*exp(-paramhat(2)*(paramhat(4)));
    else
        % Param 4 is y
        pm(5) = paramhat(1)*exp(-paramhat(2)*(paramhat(4) - paramhat(3)));
    end
    infR = pm(1)/pm(5);
    
    % Get time in years and reorder so towards the present 
    tpres = 1993;
    maxPast = tcoal(end);
    tz = linspace(tpres-maxPast, tpres, 1000);
    tpast = 1895;
    
    % Parameter estimates from Pybus2003
    if nold == 68
        % Data from 68 sequences
        oli.NC = [4095 10310 18960];
        oli.r = [0.075 0.264 0.620];
        oli.x = [1941 1953 1966];
        oli.y = [1924 1934 1943];
        oli.NA = [153 245 345];
        oli.TMRCA = [1258 1374 1481];
    else
        % Data from 63 sequences
        oli.NC = [3323 8779 15780];
        oli.r = [0.072 0.237 0.564];
        oli.x = [1941 1953 1966];
        oli.y = [1922 1932 1940];
        oli.NA = [99.6 170 251];
        oli.TMRCA = [1673 1710 1747];
    end
    
    % Get population optimal and bounds from last posterior
    idq = round(linspace(1, length(tz), 20));
    [Nf1, Nbnd1] = getPopForward(tz, qnlast, fn, tpres);
    
    % Population optimal with 95% bounds
    [Nf, Nbnd, qcum] = getPopForwardQuan(tz, qnlast, fn, tpres);
    
    % Get TMRCA from tree back in years
    TMRCA = round(tpres - maxPast);
    disp(['TMRCA = ' num2str(TMRCA)]);
    
    % PAT treatment period
    patx = [1920 1920 1980 1980];
    Nmin = min(min(Nbnd));
    Nmax = max(max(Nbnd));
    paty = [Nmin Nmax Nmax Nmin];
    
    % Get the individual parameters in this form from marginal distribs
    est.NC = [xstdlb(end, 1), xhatmean(end, 1), xstdub(end, 1)];
    est.NC = round(est.NC);
    est.r = [xstdlb(end, 2), xhatmean(end, 2), xstdub(end, 2)];
    est.x = tpres - [xstdub(end, 3), xhatmean(end, 3), xstdlb(end, 3)];
    est.x = round(est.x);
    if fnid == 5
        est.y = tpres - [xstdub(end, 4), xhatmean(end, 4), xstdlb(end, 4)];  
    else
        % Used deviation on y as deviation on y-x + mean x (approx)
        est.y = xhatmean(end, 3) + [xstdub(end, 4), xhatmean(end, 4), xstdlb(end, 4)];
        est.y = tpres - est.y;
    end
    est.y = round(est.y);
    
    % NA estimates more complicated only give mean
    if fnid == 5
        NAmean = est.NC(2)*exp(-est.r(2)*(est.x(2) - est.y(2)));
    else
        NAmean = est.NC(2)*exp(-est.r(2)*(xhatmean(end, 4)));
    end
    est.NA = round(NAmean);
    est.TMRCA = TMRCA;
    
    % Convert coalescent times to years and get a vector to overlay on PAT
    if size(tcoal, 2) == 1
        % Set tcoal to a row vector if it isn't already
        tcoal = tcoal';
    end
    yr_coal = tpres - tcoal;
    n_coal = ones(1, n);
    vx_coal = repmat(yr_coal, 2, 1);
    v_coal = [zeros(1, n); xstdub(end, 1)*ones(1, n)];
    
    % Extract oli parameters - means
    xo = oli.x(2);
    yo = oli.y(2);
    ro = oli.r(2);
    NAo_notcalc = oli.NA(2);
    NCo = oli.NC(2);
    NAo = NCo.*exp(-ro.*(xo-yo));
    
    % Extract oli std deviations - 95% treated as 2 stds
    xo_s = 0.5*(oli.x(3) - oli.x(1));
    yo_s = 0.5*(oli.y(3) - oli.y(1));
    ro_s = 0.5*(oli.r(3) - oli.r(1));
    NCo_s = 0.5*(oli.NC(3) - oli.NC(1));
    
    % Get oli curve for plotting
    I1o = tz >= xo;
    I2o = tz > yo & tz < xo;
    I3o = tz <= yo;
    Noli = NCo.*I1o + NAo.*exp(-ro.*(yo - tz)).*I2o + NAo.*I3o;
    
end

%% Visualise data - no MSE calculation possible
    
% Plot the N estimates and coalescent rate lambda
if plotindiv
    
    % Plot the tree if ultrametric from R8
    if getTree && useR8
        tr = phytreeread(name);
        h = plot(tr);
        % Get xticks and re-evaluate from TMRCA
        xt = h.axes.XTick;
        xtnew = TMRCA + xt;
        xtnew = xtnew';
        set(gca,'XTickLabel', num2str(xtnew));
        xlabel('years forward from TMRCA');
        title(['Ultrametric tree for nSeq = ' num2str(nold)]);
    end
    
    % Plot best N (assuming posterior)
    figure;
    semilogy(t, Nf2, 'k', 'LineWidth', 2);
    hold on
    semilogy(t, Nbnd2, 'g--', 'LineWidth', 1);
    semilogy(tn, Nhat, 'ro');
    hold off
    legend('final post', 'final bnds', 'iterative', 'location', 'best');
    xlabel('time');
    ylabel(['N(t) = ' fn.name]);
    title(['N(t) final estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
    xlim([tnmin tnmax]);
    
    % The parameters estimated and the std deviation bounds
    figure;
    for i = 1:numRV
        subplot(numRV, 1, i);
        plot(tn, xhatmean(:, i), 'k', tn, xstdlb(:, i),...
            'b', tn, xstdub(:, i), 'b');
        xlabel('time');
        %ylabel(['x_' num2str(i) ' in N(t) =' fn.name]);
        ylabel(['x_' num2str(i)]);
        title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        legend('estimate', '-2 std', '+2 std', 'location', 'best');
        xlim([tnmin tnmax]);
    end
    
    % Posteriors and final estimates
    figure;
    for i = 1:numRV
        subplot(numRV, 1, i);
        plot(xset{i}, qnp{i}(1, :), 'b-', xset{i}, qnp{i}(end, :), 'ko-',...
            xset{i}, qnp{i}(2:end-1, :), 'g--');
        hold on
        plot([paramhat(i) paramhat(i)], [0 max(max(qnp{i}))], 'k');
        hold off
        % Condition to account for variables with only one value
        if min(xset{i}) ~= max(xset{i})
            xlim([min(xset{i}) max(xset{i})]);
        end
        xlabel(['x_' num2str(i)]);
        ylabel(['P(x_' num2str(i) '|data)']);
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: param est = ' num2str(paramhat(i))]);
    end
    
    % Pybus2003 plot type
    if fnid == 5 || 6
        figure;
        fill(patx, paty, 'c:');
        hold on
        plot(tz, Nf, 'k', 'LineWidth', 2);
        plot(tz, Nbnd, 'r', 'LineWidth', 2);
        plot(vx_coal, v_coal, 'k:');
        hold off
        legend('PAT period', 'N(t)', 'lower bound', 'upper bound', 'location', 'best');
        xlim([tpast-100 tpres]);
        xlabel('year');
        ylabel('estimated effective number of infections');
        title('Estimated demographics of Egyptian HCV');
        
        % Version for poster
        figure;
        fill(patx, paty, 'c:');
        hold on
        plot(tz, Nf, 'k', 'LineWidth', 2);
        plot(tz, Nbnd, 'r', 'LineWidth', 2);
        %plot(vx_coal, v_coal, 'k:');
        %plot(tz, Noli, 'm', 'LineWidth', 2);
        hold off
        legend('PAT period', 'N(t)', 'lower bound', 'upper bound', 'location', 'best');
        xlim([tpres-100 tpres]);
        xlabel('year');
        ylabel('estimated effective number of infections');
        title('Estimated demographics of Egyptian HCV');
        grid
        
        % Another visualisation of the bound in relative terms
        NN = Nbnd./repmat(Nf, [2 1]);
        figure;
        plot(tz, NN);
        xlim([tpast-100 tpres]);
        xlabel('year');
        ylabel('bound(N(t))/N(t)');
        title('Relative confidence in N(t) for Egyptian HCV');
        
        % Version for poster of posterior - assume fnid == 6
        if fnid == 6
            % Pybus mean values
            olipm(1) = oli.NC(2);
            olipm(2) = oli.r(2);
            olipm(3) = oli.x(2);
            olipm(4) = oli.x(2) - oli.y(2);
            olipm = olipm';
            
            % Pybus ranges
            ranOli{1} = kron([oli.NC(1) oli.NC(3)], [1 1]);
            ranOli{2} = kron([oli.r(1) oli.r(3)], [1 1]);
            ranOli{3} = kron([oli.x(1) oli.x(3)], [1 1]);
            ranOli{4} = oli.x(2) - kron([oli.y(3) oli.y(1)], [1 1]);
            
            figure;
            % Set x in years
            zset = xset;
            zset{3} = tpres - zset{3};
            paramhatz = paramhat;
            paramhatz(3) = tpres - paramhat(3);
            for i = 1:numRV
                subplot(ceil(numRV/2), 2, i);
                hold on 
                fill(ranOli{i}, [0 max(max(qnp{i})) max(max(qnp{i})) 0], [1,1,0.9]);
                %plot(zset{i}, qnp{i}(1, :), 'b-', zset{i}, qnp{i}(end, :), 'ko-',...
                    %zset{i}, qnp{i}(2:end-1, :), 'g--');
                plot(zset{i}, qnp{i}(1, :), 'b-', 'Linewidth', 2);
                plot(zset{i}, qnp{i}(end, :), 'rs-', 'Linewidth', 2);
                plot([paramhatz(i) paramhatz(i)], [0 max(max(qnp{i}))], 'r', 'LineWidth', 2);
                plot([olipm(i) olipm(i)], [0 max(max(qnp{i}))], 'k', 'LineWidth', 2);
                hold off
                % Condition to account for variables with only one value
                if min(zset{i}) ~= max(zset{i})
                    xlim([min(zset{i}) max(zset{i})]);
                end
                xlabel(['x_' num2str(i)]);
                ylabel(['P(x_' num2str(i) '|data)']);
                legend('pybus range', 'prior', 'posterior', 'snyder estimate', 'pybus estimate', 'location', 'best');
                title(['x_' num2str(i) ': [snyder pybus] = [' num2str(round(paramhatz(i), 4, 'significant')) ', ' num2str(olipm(i)) ']']);
                grid
            end
        end
    end
end

% Clock time set to minutes
tIter = toc;
tIter = tIter/60;
disp(['Simulation time is ' num2str(tIter(end))]);
disp('********************************************************************');
for i = 1:numRV
    disp(['Estimate x_' num2str(i) ' = ' num2str(paramhat(i))]);
end
if fnid == 5 || 6
    % Additional dependent parameter calculated and ratio of infections
    disp(['Estimate x_5 = ' num2str(pm(5))]);
    disp(['Infection ratio = ' num2str(infR)]);
    if fnid == 5
        disp(['Length of exponential growth = ' num2str(pm(4) - pm(3))]);
    else
        disp(['Length of exponential growth = ' num2str(pm(4))]);
    end
end
    
% Store individual run removing heavy, unneeded variables
save(['tree63True_' num2str(numRV)], 'fn', 'nLin', 'tn', 'lamhat', 'Nhat', 'statsStruc',...
    'tIter', 'tcoal', 'name', 'remDuplicates', 'treeSeq', 'evDist');