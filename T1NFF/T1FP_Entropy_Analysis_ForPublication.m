clear all;
close all;
addpath('../function/');
base_dir = uigetdir;

%% For All chronic timepoints
% Weighted, 40 bins
T1_entropy = [3.292015467	3.783864573	2.973529192	2.563597534	3.623484896	4.039070068	3.521463398	3.434995204	3.938055327	3.762232561	3.277875859	3.242358077	2.960718726	2.971889596	2.802255103	2.869634886	3.151586898	3.086717588	3.058121378	3.074957238	2.578228703	3.295680752	3.382698877	2.956554086	2.763665867	2.595815638	3.213119346	2.910791296	3.027350691	3.415119561	3.049826019	2.827704325	3.329855847	2.252669124	2.427804709	2.078900194	1.785305929	1.986330983	2.162572264	2.073256383	2.641822704	2.381059023	2.243795798	2.961622881	3.009159374	2.27613887	1.481916078	2.855204435	2.767544556	2.197571944	2.328955153	3.447057878	3.819486644	3.364619	3.415637825	3.653883869	3.782141879	3.519648062	3.260072675	3.691184573	3.545920887	3.352747155	2.638744434	3.081907058	2.755593336	2.521004546	1.905956725	2.476013596	2.486637372	2.574397452	2.262907233	3.353903554	3.636646253	2.874416194	3.267397732	3.408997009	2.918151106	3.577349186	3.83306413	3.491983393	2.594694996	3.672898514	3.262030675	3.193682288	2.792275342	3.491605214	3.137371244	3.330624871	2.71342767	3.590609451	3.364294755	3.349767401	3.569175429	3.607389089	2.97435778	2.321653449	2.419381946	2.860320778];


% Remove 2 pixels
% T1_entropy = [3.226650665	3.457850828	2.718079221	2.499056589	3.371479193	3.195670399	3.064785548	2.804755242	2.975664237	3.392783749	3.377906439	3.008738456	2.98838875	2.844429901	2.832724493	2.697805734	2.822762112	3.029865506	2.81136788	2.954446942	2.725400327	2.409891914	2.474903441	2.778625884	2.656873042	2.252382635	2.70835284	3.030668964	2.639595553	2.359077838	3.071203999	2.84446803	2.436041881	2.763742186	1.992180025	2.203703282	1.993940438	1.785305929	1.986330983	2.02325357	1.998485279	2.457795368	2.290983872	2.198052043	2.581775734	2.678328269	2.115803637	1.332820405	2.516182382	2.625713045	2.010239527	1.914737636	3.197997566	3.505325142	3.074497691	3.089220972	2.908192061	2.84738103	2.937687855	2.596840402	3.000211946	3.169069808	2.994413901	2.524967107	2.998447499	2.565253142	2.20212135	1.835902812	2.329831286	2.427465917	2.198021807	1.87599672	3.175221734	3.322749077	2.377708927	3.182048804	2.95398955	2.545954336	2.951999651	3.290760029	3.080862045	2.380860235	3.293971748	3.007071329	2.189238256	2.485228136	3.241956654	2.850966975	2.509007051	1.557431962	3.330939307	3.090225788	3.0843122	3.30905914	3.409117229	2.773133387	1.959298911	0	1.554585169];

% Remove 1 pixel
% T1_entropy = [3.292015467	3.521195067	2.79463277	2.499056589	3.445044499	3.621144808	3.297626439	2.892417091	3.066471628	3.462203302	3.592682664	3.177138664	3.0530908	2.897942434	2.896518766	2.697805734	2.822762112	3.085583319	2.92440905	3.019313842	2.947289192	2.487492159	2.979806954	3.133844038	2.753393644	2.505249493	2.474903441	2.778625884	2.812441399	2.815071955	3.284571717	3.049826019	2.737397331	3.024822002	1.992180025	2.203703282	1.993940438	1.785305929	1.986330983	2.111291841	1.998485279	2.572645343	2.290983872	2.198052043	2.758429824	2.932287882	2.115803637	1.332820405	2.516182382	2.691009153	2.010239527	2.105397711	3.353331226	3.701725673	3.228027809	3.203343525	3.465365569	3.524141992	3.210316365	3.018648943	3.405099327	3.321758834	3.099583357	2.59655574	2.998447499	2.652450464	2.373915237	1.835902812	2.329831286	2.427465917	2.436638708	1.87599672	3.175221734	3.41357256	2.694810802	3.267397732	3.105601702	2.641704474	3.344269698	3.454493031	3.199919752	2.380860235	3.376717816	3.097800573	2.425098369	2.672905595	3.373857357	2.953169881	2.980240673	2.127065631	3.449855635	3.235358758	3.349767401	3.374464205	3.607389089	2.928387789	1.959298911	0.970950594	1.924174352];

% unweighted, 40 bins
% T1_entropy = [3.537273219	4.032394521	3.326644545	3.008024795	3.927864033	4.258540624	3.857469858	3.479601779	3.799957565	4.233525404	4.002220303	3.694694354	3.571369543	3.598695476	3.482047247	3.367627438	3.436607051	3.43726271	3.513891595	3.47255582	3.540584069	3.085670113	3.288293569	3.593676212	3.259061911	2.704099781	2.617536533	3.276420835	3.067476457	2.652345561	3.928134475	3.478835254	3.43127238	3.768273847	2.4758585	2.74216951	2.303627993	2.232848992	2.383988609	2.682792823	2.516696304	3.079797834	2.906915359	2.826764657	3.365261528	3.604876133	2.936911074	2.214199695	3.118697829	3.342514592	2.697405766	2.556359143	3.502596841	4.049195895	3.654991497	3.415678079	3.694891829	4.018089708	3.975646384	3.30708252	4.047772969	3.940762991	3.723397038	3.237951205	3.56244661	3.198341371	2.739881706	2.271066057	2.561964212	2.893388358	2.879451104	2.765574736	3.190119875	3.52847719	3.197407646	3.392662215	3.714523193	3.109386672	3.564926018	4.125803243	3.501826416	2.928782852	4.001909518	3.65633558	3.374548236	2.96147252	3.995461119	3.552783151	3.358922883	3.269576746	3.66462172	3.290420494	3.343119911	3.55559154	3.427307576	3.022182094	2.307820398	2.05881389	3.090032777];

% unweighted remove 1 pixels
% T1_entropy = [3.537273219	3.782380131	3.286405173	2.888662706	3.84755692	3.864413376	3.777840376	3.432148197	3.366493661	3.944728716	3.725613275	3.458720552	3.394581786	3.484402087	3.447749151	3.226404779	3.352014857	3.310557034	3.394060013	3.436646427	3.463385702	3.044935395	3.137126814	3.294681477	2.934717607	2.532314007	2.36559623	3.029408653	2.772208352	2.41268789	3.692190757	3.421283068	3.356971263	3.56895206	2.288466328	2.529550926	2.220323638	2.15218353	2.383988609	2.588147006	2.446481784	3.014843837	2.781837408	2.785245314	3.273114274	3.604876133	2.796907293	2.086641775	2.962271231	3.308123598	2.697405766	2.236779563	3.3193688	3.898256069	3.570623578	3.131673258	3.239668316	3.617076226	3.840960127	3.074820572	3.896728115	3.73576803	3.494882387	3.08415993	3.486437276	3.103631795	2.64618737	2.271066057	2.354435538	2.83859057	2.751154309	2.560352492	3.065543713	3.360131473	3.030869578	3.340515087	3.388982986	2.898502483	3.324424043	3.818738355	3.140864629	2.627413428	3.824450462	3.554838897	2.503997527	2.727827978	3.89027626	3.327108655	3.103347507	2.935649892	3.525242792	3.247439407	3.303235567	3.479055348	3.427307576	2.975447368	1.94400975	0.918295834	2.235926351];

% Unweighted, 20 bins
T1_entropy = [2.574698565	3.11641416	2.358427424	2.111429903	2.978402775	3.347837177	2.96493762	2.894477301	3.259330177	3.115440133	2.766654228	2.611181993	2.602057079	2.524576516	2.406260566	2.496369982	2.538348015	2.509035772	2.48490349	2.585283285	2.148374255	2.398838914	2.648555312	2.402159542	1.839126497	1.622556249	2.368758408	2.256482054	1.903881788	2.997549072	2.689013128	2.520381504	2.902648639	1.608617724	1.961731166	1.490580398	1.360946274	1.447984595	1.83254541	1.572194612	2.203974463	1.919007873	1.870708629	2.587282248	2.816679246	2.041220906	1.243180398	2.188882993	2.38232161	1.800000911	1.729393449	2.549083431	3.112251948	2.690413038	2.448794937	2.799663983	3.151489355	3.076284397	2.523777569	3.121384648	3.038978984	2.781351419	2.325361734	2.543999705	2.297340406	1.881329453	1.220375637	1.62548643	1.93954513	2.001699743	1.835828979	2.509446836	2.606201621	2.257344893	2.549007671	2.747476322	2.19547983	2.732843101	3.221968463	2.621932649	2.187450335	3.111954791	2.733228663	2.609931016	2.162927762	3.004677481	2.661841933	2.472334338	2.31850778	2.737490158	2.371470134	2.427307582	2.619791797	2.490255144	2.080386854	1.717040781	1.224394445	2.294274067];
T1_entropy = [2.993392576	3.492984012	2.862405476	2.457504452	3.39451505	3.740589706	3.346396138	3.335937138	3.72325104	3.562172475	3.157377062	3.129883467	3.05689881	3.005598711	2.882839099	2.988001102	2.941133856	2.997195895	2.971977065	3.055517486	2.626664186	2.835322898	3.085450382	2.735405555	2.250710988	2.173542317	2.71618431	2.677010072	2.324117032	3.491558497	3.104346786	2.801969921	3.294875277	2.044585772	2.329539081	1.742622512	1.869770343	1.977758629	2.186478613	2.022071551	2.554632948	2.385828165	2.31047321	2.962999385	3.186802793	2.421355088	1.540923456	2.616716926	2.88429002	2.237278746	2.103219547	3.059415774	3.537625908	3.163406315	2.921730336	3.219368163	3.499090921	3.485181749	2.913227995	3.533257224	3.440154081	3.248961176	2.747676423	3.030186078	2.686528728	2.127363525	1.719653421	2.035525092	2.423369853	2.402595621	2.222124456	2.88896271	3.026675242	2.620726767	2.896002098	3.201135741	2.619286961	3.128480775	3.606875407	3.023956088	2.640707151	3.470379829	3.20414875	3.010428829	2.567099536	3.434291389	3.191528775	2.97644311	2.757543733	3.208378232	2.738875143	2.881210404	3.000204708	2.930433724	2.49402299	1.962354484	1.752715279	2.610063541	2.094856784	1.162563308];

T1_entropy_remote = [0.865856617	0.67694187	1.480290753	1.2586633	1.236386411	1.423794941	1.681046078	1.251629167	1.494918848	1.337778097	1.700438513	1.323489961	1.58468025	1.672492724	1.429152984	2.134711144	1.952819531	1.334599426	1.620538846	1.901872151	1.54536616	1.846393983	1.987773371	1.945661003	2.022887961	1.294995075	1.422475118	1.573114618	1.916445381	1.040852083	1.481863125	1.696842803	1.870018047	0.985228136	0.99631652	1.253297578	2.230337797	1.530639062	1.086312841	0.976020648	1.392147224	2.43147771	1.726474118	1.231724487	0.811847521	1.485475297	1.021119189	0.811278124	1.337306179	0.73206669	1.298794941	1.052981575	0.937301241	0.964086328	1.303169878	1.224394445	0.787004656	1.358584608	1.632105718	1.376074607	1.924820746	2.095646899	1.217016339	0.793992434	1.030426461	0.970950594	1.56909404	0.959686894	1.722971279	1.226817429	1.79650626	1.432928309	2.09790657	1.571091977	1.742966195	1.646780976	1.67673703	1.214540408	1.468590448	1.254999238	1.130329644	1.882045108	1.217038755	1.109227139	1.115644407	1.469182539	1.570950594	1.421554169	1.248609233	2.087210311	1.781839569	1.906305535	0.991264261	1.514723983	0.959686894	1.325784781	1.498395816	1.467457965	0.266764988	1.582683189];

FF_mean = [8.67898015	6.254385141	3.163614774	2.668638449	10.8550361	14.04099999	16.97053298	15.60447728	14.78351024	12.17616673	8.138109019	6.817493369	4.872403662	3.686957353	3.121698421	3.789357016	6.199221444	6.298566453	6.349978604	3.980585046	4.247593918	4.548496437	5.296721288	5.879011562	1.701985082	4.852949747	4.657240867	5.177199829	1.631045615	1.617945671	2.868202953	0.57090047	2.862495792	2.728173187	2.586949272	2.267199355	1.965339977	2.662233876	3.587132183	3.773357303	3.1280887	3.132784497	3.067722853	3.436115497	3.447189287	2.551688089	1.972721009	4.457039614	3.453607224	1.485206799	9.884391748	6.97808652	5.563248478	6.763654469	9.737753159	10.52630823	9.897422831	10.19417579	10.56805944	11.57662087	12.4761005	6.204056188	6.20985808	7.109428801	8.894293271	3.431069405	2.016370513	1.012086171	5.763136511	8.114717589	4.279872056	5.791892352	6.859122164	4.279127297	7.951895703	2.760749916	3.46983488	8.439058352	12.71780664	3.220006784	2.680067088	10.1953705	7.028997074	4.972657679	4.079221614	9.364701114	7.498762478	3.849234679	3.67828082	3.011781264	3.011781264	1.87500241	2.795359749	3.624650895	2.636651685	1.430597131	3.262081827	4.778726584];
FF_mean = [8.67898015	6.254385141	3.163614774	2.668638449	10.8550361	14.04099999	16.97053298	15.60447728	14.78351024	12.17616673	8.138109019	6.817493369	4.872403662	3.686957353	3.121698421	3.789357016	6.199221444	6.298566453	6.349978604	3.980585046	4.247593918	4.548496437	5.296721288	5.879011562	1.701985082	4.852949747	4.657240867	5.177199829	1.631045615	1.617945671	2.868202953	0.57090047	2.862495792	2.728173187	2.586949272	2.267199355	1.965339977	2.662233876	3.587132183	3.773357303	3.1280887	3.132784497	3.067722853	3.436115497	3.447189287	2.551688089	1.972721009	4.457039614	3.453607224	1.485206799	9.884391748	6.97808652	5.563248478	6.763654469	9.737753159	10.52630823	9.897422831	10.19417579	10.56805944	11.57662087	12.4761005	6.204056188	6.20985808	7.109428801	8.894293271	3.431069405	2.016370513	1.012086171	5.763136511	8.114717589	4.279872056	5.791892352	6.859122164	4.279127297	7.951895703	2.760749916	3.46983488	8.439058352	12.71780664	3.220006784	2.680067088	10.1953705	7.028997074	4.972657679	4.079221614	9.364701114	7.498762478	3.849234679	3.67828082	3.011781264	3.011781264	1.87500241	2.795359749	3.624650895	2.636651685	1.430597131	3.262081827	4.778726584 2.18 1.30];

r2star_mean = [61.66396949	66.72801911	50.38543069	49.13664347	61.92823034	67.09022474	56.0265886	62.97160686	61.44405532	54.40868886	48.09030404	41.94485716	36.68863286	61.51613479	34.331597	35.26857996	49.84719683	44.13125753	48.86915895	33.58955265	38.02172821	43.66634368	45.68736055	51.37058864	41.99326427	41.81526572	39.76389763	44.32598027	40.20305723	20.62063077	34.15637537	32.2809834	33.40768945	36.12753645	33.61147748	36.66169398	25.66525328	38.247155	33.68132133	35.85678127	37.84415544	35.12353681	39.35319758	37.00863445	36.9487465	36.67803083	37.61043845	45.06678416	37.977241	33.36696338	63.3389743	65.71483958	72.28083434	56.24927084	54.09668959	51.64011137	51.37432136	46.7566622	55.72579252	50.11184238	56.23618837	51.96202513	45.32157871	53.42367437	46.13076239	44.17880923	43.03302699	35.28729478	38.3749292	43.1041654	36.08635756	47.53581232	48.63565039	46.02452838	38.34390915	40.34869025	39.22262539	104.1659527	89.32698057	68.73375711	46.10527221	66.09715616	59.76106088	50.18403992	45.89983436	66.9378156	56.83786406	48.12828032	35.04357941	48.68408956	48.68408956	43.02469319	50.1362844	47.87491294	43.30856064	32.83326505	32.14249061	35.03343475];
r2star_mean = [61.66396949	66.72801911	50.38543069	49.13664347	61.92823034	67.09022474	56.0265886	62.97160686	61.44405532	54.40868886	48.09030404	41.94485716	36.68863286	61.51613479	34.331597	35.26857996	49.84719683	44.13125753	48.86915895	33.58955265	38.02172821	43.66634368	45.68736055	51.37058864	41.99326427	41.81526572	39.76389763	44.32598027	40.20305723	20.62063077	34.15637537	32.2809834	33.40768945	36.12753645	33.61147748	36.66169398	25.66525328	38.247155	33.68132133	35.85678127	37.84415544	35.12353681	39.35319758	37.00863445	36.9487465	36.67803083	37.61043845	45.06678416	37.977241	33.36696338	63.3389743	65.71483958	72.28083434	56.24927084	54.09668959	51.64011137	51.37432136	46.7566622	55.72579252	50.11184238	56.23618837	51.96202513	45.32157871	53.42367437	46.13076239	44.17880923	43.03302699	35.28729478	38.3749292	43.1041654	36.08635756	47.53581232	48.63565039	46.02452838	38.34390915	40.34869025	39.22262539	104.1659527	89.32698057	68.73375711	46.10527221	66.09715616	59.76106088	50.18403992	45.89983436	66.9378156	56.83786406	48.12828032	35.04357941	48.68408956	48.68408956	43.02469319	50.1362844	47.87491294	43.30856064	32.83326505	32.14249061	35.03343475 28.95880155 30.29987874];

%LGE_entropy = [3.919292542	4.110784112	4.089974654	3.810200708	3.846419175	3.924069395	3.902059585	3.870858235	4.07271552	4.109600022	3.991972408	3.55626904	3.716552049	3.58165343	3.801258766	3.861091469	3.903447958	3.857795248	3.924375674	3.998552725	3.854297105	3.630008655	3.716048169	3.787515764	3.744452936	3.800961858	3.708158957	4.109147895	3.912950805	3.658307272	3.579709702	3.411269011	3.524858076	3.909323454	3.976122474	3.766064815	3.732822447	3.936883259	3.714602979	3.856571364	3.727014989	3.920935007	3.500835396	3.921629125	3.744673551	3.868735721	3.656300872	3.818267946	3.894800055	3.535505743	4.019492796	3.862468727	3.89402203	3.812185276	3.744416443	3.893930839	3.992244852	3.918095732	3.866078746	3.97073485	3.944800313	3.789156797	3.797879395	3.687910764	3.712346131	3.843611986	3.97440545	3.669000969	3.722775933	4.041433827	3.949428307	3.929890893	3.889734359	3.837547193	3.835833493	4.001945187	4.039818756	3.485067125	3.972056199	4.031901369	4.016952846	3.504339045	3.984042202	3.921238236	3.834988021	3.093944588	3.769546283	4.089783055	3.858992869	4.149593829	4.05998882	4.010941276	4.024992838	4.024992838	3.924686676	3.771527784	3.236857303	3.488525295];
%T2_mean = [35.73141593	33.65185185	36.2531401	37.49791667	33.84610169	31.33308824	33.07565982	34.95344828	30.67777778	30.97619048	32.88369565	37.26329114	36.15433526	39.67062147	41.66467066	35.97469136	37.01654676	35.832	36.44820144	40.23121693	38.1	36.75977011	34.11805556	36.82105263	39.25058824	30.67222222	35.85063291	35.27901235	32.83033708	37.63902439	34.94651163	34.87884615	42.80519481	35.74	37.31527778	36.992	38.79032258	41.92209302	42.43303965	39.26888889	39.3828125	38.70470085	38.43084112	36.57607362	34.93373494	37.15578231	41.98053691	45.42808989	50.81202346	44.57482517	34.77634409	34.2490099	31.89351351	31.71322751	33.26290323	33.10379747	29.2627907	32.53714286	33.30759494	30.80294118	32.09134615	33.48412698	38.36223776	35.68367347	37.07824561	34.49954751	31.49328859	29.97372263	34.46447368	37.78611111	39.50854701	37.10631579	36.1956044	37.21166667	34.87608696	34.23253968	44.94833333	39.47894737	33.06161616	38.12216216	38.5125	37.14335664	33.948	38.85479452	34.64666667	44.73380952	33.00943396	38.7976	33.37536232	35.14545455	34.14265734	32.81756757	32.47	32.228125	32.92677165	33.61666667	37.98518519	39.2];
%T1_mean = [1188.030534	1231.571429	1244.289063	1253	1238.240741	1210.397849	1239.552632	1247.231579	1240.596774	1232.391304	1272.680412	1230.860927	1349.233503	1317.391026	1315.461538	1340.615385	1234.606557	1266.746032	1215.628378	1269.82963	1238.860294	1157.254545	1191.208333	1201.549296	1136.527273	1171.294118	1152.04	1173.5	1128.472222	1398.142857	1344.540541	1357.68	1304.037037	1132.917526	1161.135135	1121.179104	1137.507246	1236.233645	1223.196721	1214.12987	1224.135135	1217.198582	1233.151079	1257.588785	1262.894737	1357.328358	1256.8	1202.6	1312.226415	1230.882353	1184.111111	1171.688679	1221.016667	1218.283186	1114.542373	1132.634921	1202.987805	1278.043478	1138.113208	1217.135135	1212.564815	1196.351351	1267.375	1297.909774	1286.614679	1210.398438	1199.294737	1173.230769	1236.350515	1224.948718	1197.864865	1127.302632	1143.594937	1164.363636	1172.08046	1176.710744	1177.376238	1194.680851	1218.494737	1224.775862	1366.738095	1290.521277	1295.114943	1214.433333	1246.419355	1365.908397	1325.428571	1197.956522	1228.566667	1125.738255	1141.071429	1166.573529	1158.844444	1127.158333	1159.008475	1147.136364	1141.7	1248.105263];
dfp_treated_label = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1,...
    1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0];
hemo_label = [1	1	1	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	1	1	0	0	1	0	1	1	1	1	1	1	1	1	1	0	1	1	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	1 0 0];

r2star_mean_core = [61.94307828	76.18135322	59.87638549	0	64.67146453	68.18503134	62.58521901	0	63.074596	69.00121096	58.57586781	52.44060189	45.65815824	41.06407289	76.54461979	40.30648736	42.38267662	46.64519181	56.46477305	54.71190635	38.51559833	39.23125351	44.43670648	48.21565244	46.20186406	43.87674857	47.8662608	53.48974484	55.13819647	44.42573949	0	0	35.08285907	57.72704076	0	0	51.00690437	0	41.55915137	32.26982206	40.97986235	40.49371765	33.70151854	47.53551184	49.40957938	41.55422674	41.20528461	0	59.29470013	47.65263222	0	66.70692886	74.22914318	76.70050628	68.2719432	72.67829265	55.37004556	59.03588732	49.5543491	61.6618463	57.35343764	61.51215175	60.57790793	52.59371398	63.96942263	70.62443949	49.9444309	44.72319865	60.33931248	41.21062379	50.15432889	43.08639105	58.80959883	59.45235923	53.23263297	49.41202765	47.05307283	43.88494615	104.1659527	105.3735606	71.3822877	55.18135139	87.10114117	77.5029901	58.72351064	50.30429298	71.734744	58.91634099	55.35176295	37.33161715	62.2854226	0	55.91131405	56.37855216	57.49797238	52.77414554	0	0	38.99169973 0 0];

name_label = {'Merry_6MO', 'Merry_6MO', 'Merry_6MO', 'Merry_6MO', 'Merry_1YR', 'Merry_1YR', 'Merry_1YR', 'Merry_15YR', 'Merry_15YR', 'Merry_15YR', 'Merry_15YR',...
    'Ryn_9MO', 'Ryn_9MO', 'Ryn_9MO', 'Ryn_9MO', 'Ryn_9MO', 'Ryn_1YR', 'Ryn_1YR', 'Ryn_1YR', 'Ryn_1YR', 'Ryn_1YR',...
    'Mojave_9MO', 'Mojave_9MO', 'Mojave_9MO', 'Mojave_9MO', 'Mojave_1YR', 'Mojave_1YR', 'Mojave_1YR', 'Mojave_1YR',...
    'Sahara_9MO', 'Sahara_9MO', 'Sahara_9MO', 'Sahara_9MO', 'Sahara_1YR', 'Sahara_1YR', 'Sahara_1YR', 'Sahara_1YR',...
    'ZZ_6MO', 'ZZ_6MO', 'ZZ_6MO', 'ZZ_6MO', 'ZZ_1YR', 'ZZ_1YR', 'ZZ_1YR', 'ZZ_1YR',...
    'Tina_6MO', 'Tina_6MO', 'Tina_1YR', 'Tina_1YR', 'Tina_1YR',...
    'Sunny_6MO', 'Sunny_6MO', 'Sunny_6MO', 'Sunny_6MO', 'Sunny_1YR', 'Sunny_1YR', 'Sunny_1YR', 'Sunny_1YR', 'Sunny_15YR', 'Sunny_15YR', 'Sunny_15YR', 'Sunny_15YR'};

find(abs(zscore(T1_entropy))>3)
find(abs(zscore(FF_mean))>3)     % 7, 8
find(abs(zscore(r2star_mean))>3) % 79, 80
% find(abs(zscore(LGE_entropy))>3) % 87, 98
% find(abs(zscore(T2_mean))>3) % 50

find(abs(zscore(T1_entropy))>2.5)  % 48
find(abs(zscore(FF_mean))>2.5)     % 7, 8, 9
find(abs(zscore(r2star_mean))>2.5) % 79, 80
% find(abs(zscore(LGE_entropy))>2.5) % 87, 98
% find(abs(zscore(T2_mean))>2.5) % 50

find(abs(zscore(T1_entropy))>2)  % 38, 48, 68
find(abs(zscore(FF_mean))>2)     % 6, 7, 8, 9, 10
find(abs(zscore(r2star_mean))>2) % 31, 79, 80
% find(abs(zscore(LGE_entropy))>2) % 33, 87, 98
% find(abs(zscore(T2_mean))>2) % 49, 50, 51, 78, 87

exclude_idx_ff = zeros(1, length(FF_mean)) > 0;
exclude_idx_r2star = zeros(1, length(r2star_mean)) > 0;
exclude_idx = zeros(1, length(FF_mean)) > 0;


% exclude_idx_ff = FF_mean < 4;
% exclude_idx_r2star = r2star_mean < 35; % T1_entropy = 2.5968
% 
exclude_idx_ff = FF_mean >= 4;
exclude_idx_r2star = r2star_mean >= 35; % T1_entropy = 2.1776

% exclude_idx_ff = FF_mean < 4;
% exclude_idx_r2star = r2star_mean >= 35;
% 
% exclude_idx_ff = FF_mean >= 4;
% exclude_idx_r2star = r2star_mean < 35; % T1_entropy = 2.1552

% exclude_idx_r2star = r2star_mean < 40;
exclude_idx_hemo = hemo_label;
exclude_idx = exclude_idx_ff | exclude_idx_r2star;

%exclude_idx = exclude_idx_ff | exclude_idx_r2star | exclude_idx_hemo;

% 2.5SD
% exclude_idx(7) = 1;
% exclude_idx(8) = 1;
% exclude_idx(9) = 1;
% exclude_idx(48) = 1;
% exclude_idx(87) = 1;
% exclude_idx(98) = 1;
% exclude_idx(79) = 1;
% exclude_idx(80) = 1;
% 

% % excluded by T1-entropy and FF
% exclude_idx(8) = 1;
% exclude_idx(52) = 1;
% exclude_idx(66) = 1;
% exclude_idx(71) = 1;
% 
% % excluded by T1-entropy and r2star
% exclude_idx(79) = 1;
% exclude_idx(68) = 1;
% 
% % excluded by image quality
% exclude_idx(80) = 1;
% 
% % exclude_idx(7) = 1;
% % exclude_idx(9) = 1;
% % exclude_idx(54) = 1;
% % exclude_idx(64) = 1;
% % exclude_idx(74) = 1;
% % exclude_idx(79) = 1;
% % exclude_idx(80) = 1;
% 
% 
% 
% % 
% % % excluded by T1-entropy and LGE-entropy
% % exclude_idx(87) = 1;
% % exclude_idx(98) = 1;
% % 
% % % % exclude by T1-entropy and T2-mean
% exclude_idx(50) = 1;
% % exclude_idx(49) = 1;
% % exclude_idx(51) = 1;
% % %exclude_idx(27) = 1;
% % %exclude_idx(69) = 1;
% % exclude_idx(78) = 1;
% % 
% % % exclude for T1-entropy and R2star
% exclude_idx(31) = 1;


% excluded by T1-entropy and FF
%exclude_idx(8) = 1;
exclude_idx(52) = 1;
exclude_idx(66) = 1;
exclude_idx(71) = 1;

% excluded by T1-entropy and r2star
exclude_idx(79) = 1;
exclude_idx(68) = 1;

% excluded by image quality
exclude_idx(80) = 1;

%exclude_idx(7) = 1;
%exclude_idx(9) = 1;
exclude_idx(54) = 1;
exclude_idx(78) = 1;

% After 20 bins
exclude_idx(51) = 1;
% exclude_idx(79) = 1;
% exclude_idx(68) = 1;

% exclude for FF and R2star
exclude_idx(53) = 1;

% exclude_idx(26) = 1;
% For image quality
exclude_idx(30) = 1; % Sahara 9MO slice1




T1_entropy(exclude_idx) = [];
r2star_mean(exclude_idx) = [];
FF_mean(exclude_idx) = [];
T1_entropy_remote(exclude_idx) = [];
% LGE_entropy(exclude_idx) = [];
% T2_mean(exclude_idx) = [];
% T1_mean(exclude_idx) = [];

dfp_treated_label(exclude_idx) = [];

r2star_mean_core(exclude_idx) = [];

%% Main Analysis Here
color_cell1 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell2 = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};


addpath('../function/BlandAltman/');
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'chronic'};

% figure();
% plot(T1_entropy, FF_mean, 'o');
% 
% figure(); 
% plot(T1_entropy, r2star_mean, 'o');


labels = repmat(1, [1, length(T1_entropy)]);
% Baseline data with noise
data1 = cat(3, T1_entropy(:));
data2 = cat(3, FF_mean(:));

data1 = cat(3, FF_mean(:));
data2 = cat(3, T1_entropy(:));

% BA plot paramters
tit = 'T1 Entropy vs Fat Fraction'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
% label = {'T1-Entropy','FF','%'}; % Names of data sets
label = {'FF','T1-entropy'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [0 20 0 5];
%axesLimits = [2 6 0 20];
% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'AxesLimits', axesLimits);


%%
labels = repmat(1, [1, length(T1_entropy)]);
% Baseline data with noise
data1 = cat(3, T1_entropy(:));
data2 = cat(3, r2star_mean(:));

data1 = cat(3, r2star_mean(:));
data2 = cat(3, T1_entropy(:));

% BA plot paramters
tit = 'T1 Entropy vs R2star'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'T1-Entropy','R2star','s-1'}; % Names of data sets
label = {'R2star','T1-entropy','s-1'}; % Names of data sets

corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [2 6 0 100];
axesLimits = [0 100 0 5];
% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'AxesLimits', axesLimits);

%% FF vs R2star
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'chronic'};

% figure();
% plot(T1_entropy, FF_mean, 'o');
% 
% figure(); 
% plot(T1_entropy, r2star_mean, 'o');


labels = repmat(1, [1, length(T1_entropy)]);
% Baseline data with noise

data1 = cat(3, FF_mean(:));
data2 = cat(3, r2star_mean(:));

% BA plot paramters
tit = 'Fat Fraction vs R2star'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
% label = {'T1-Entropy','FF','%'}; % Names of data sets
label = {'FF','R2star'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [0 20 0 100];
%axesLimits = [2 6 0 20];
% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'AxesLimits', axesLimits);

%% For publication
color_cell_exvivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
color_cell_invivo = {[242,240,247]/255, [203,201,226]/255, [158,154,200]/255, [117,107,177]/255, [84,39,143]/255}; % Purple
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};

FF_mean_illu = zeros(1, length(FF_mean));
FF_mean_illu_label_pp = find(FF_mean >= 4 & r2star_mean >= 35);
FF_mean_illu_label_np = find(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_label_nn = find(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu_pp = FF_mean(FF_mean >= 4 & r2star_mean >= 35) + 0.6;
FF_mean_illu_np = FF_mean(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_nn = FF_mean(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu = zeros(1, length(r2star_mean));
r2star_mean_illu_label_pp = find(FF_mean >= 4 & r2star_mean >= 35);
r2star_mean_illu_label_np = find(FF_mean < 4 & r2star_mean >= 35);
r2star_mean_illu_label_nn = find(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu_pp = r2star_mean(FF_mean >= 4 & r2star_mean >= 35) + 3;
r2star_mean_illu_np = r2star_mean(FF_mean < 4 & r2star_mean >= 35) + 3;
r2star_mean_illu_nn = r2star_mean(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu(FF_mean_illu_label_pp) = FF_mean_illu_pp;
FF_mean_illu(FF_mean_illu_label_np) = FF_mean_illu_np;
FF_mean_illu(FF_mean_illu_label_nn) = FF_mean_illu_nn;

r2star_mean_illu(r2star_mean_illu_label_pp) = r2star_mean_illu_pp;
r2star_mean_illu(r2star_mean_illu_label_np) = r2star_mean_illu_np;
r2star_mean_illu(r2star_mean_illu_label_nn) = r2star_mean_illu_nn;

figure();
scatter(FF_mean_illu, r2star_mean_illu, 'filled');
xlim([0 20]); ylim([0 100]);
hold on;
xline(4.25);
yline(36);

T1_entropy_illu_pp = T1_entropy(FF_mean >= 4 & r2star_mean >= 35);
T1_entropy_illu_np = T1_entropy(FF_mean < 4 & r2star_mean >= 35);
T1_entropy_illu_nn = T1_entropy(FF_mean < 4 & r2star_mean < 35);

% T2_mean_illu_pp = T2_mean(FF_mean >= 4 & r2star_mean >= 35);
% T2_mean_illu_np = T2_mean(FF_mean < 4 & r2star_mean >= 35);
% T2_mean_illu_nn = T2_mean(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu_pp_v2 = FF_mean(FF_mean >= 4 & r2star_mean >= 35);
FF_mean_illu_np_v2 = FF_mean(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_nn_v2 = FF_mean(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu_pp_v2 = r2star_mean(FF_mean >= 4 & r2star_mean >= 35);
r2star_mean_illu_np_v2 = r2star_mean(FF_mean < 4 & r2star_mean >= 35);
r2star_mean_illu_nn_v2 = r2star_mean(FF_mean < 4 & r2star_mean < 35);

% Exclusion
FF_mean_illu_np_v2([15 17 19]) = [];
T1_entropy_illu_np([15 17 19]) = [];
r2star_mean_illu_np_v2([15 17 19]) = [];

FF_mean_illu_np([15 17 19]) = [];
r2star_mean_illu_np([15 17 19]) = [];

%% For publication
plotHandles = zeros(4,3);
% figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(r2star_mean_illu_pp, FF_mean_illu_pp,  64, 'filled', 'Color', color_cell_exvivo{1});
hold on;
plotHandles(:,2) = scatter(r2star_mean_illu_np, FF_mean_illu_np, 64, 'filled', 'Color', color_cell_avg16{1});
plotHandles(:,3) = scatter(r2star_mean_illu_nn, FF_mean_illu_nn, 64, 'filled', 'Color', color_cell_invivo{1});

set(plotHandles(:,2), 'Marker', 's', 'Color', color_cell_exvivo{4});
set(plotHandles(:,2), 'Marker', 's',  'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_cell_exvivo{5}, 'MarkerFaceColor' , color_cell_exvivo{2});


set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

set(plotHandles(:,3), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,3), 'Marker', '^', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

ylim([0 20]); xlim([0 100]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu, FF_mean_illu, 1);
[r, p] = corrcoef(r2star_mean_illu, FF_mean_illu); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu, FF_mean_illu);
stats.corrP = p(1,2);
stats.N = length(r2star_mean_illu);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu) - FF_mean_illu).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:100;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
%axis off;



%%  Red
color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(r2star_mean_illu_pp, FF_mean_illu_pp,  96, 'filled', 'Color', color_cell_exvivo{1});
hold on;
plotHandles(:,2) = scatter(r2star_mean_illu_np, FF_mean_illu_np, 96, 'filled', 'Color', color_cell_avg16{1});
plotHandles(:,3) = scatter(r2star_mean_illu_nn, FF_mean_illu_nn, 96, 'filled', 'Color', color_cell_invivo{1});


set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

set(plotHandles(:,2), 'Marker', 's', 'Color', color_cell_gray{4});
set(plotHandles(:,2), 'Marker', 's',  'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_cell_gray{5}, 'MarkerFaceColor' , color_cell_gray{2});

set(plotHandles(:,3), 'Marker', '.', 'Color', color_cell_gray{4});
set(plotHandles(:,3), 'Marker', '^', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_gray{5}, 'MarkerFaceColor', color_cell_gray{2});

ylim([0 20]); xlim([0 100]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_pp, FF_mean_illu_pp, 1);
[r, p] = corrcoef(r2star_mean_illu_pp, FF_mean_illu_pp); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_pp, FF_mean_illu_pp);
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_pp);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_pp) - FF_mean_illu_pp).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 36:1:100;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
%axis off;

%% Blue
color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(r2star_mean_illu_pp, FF_mean_illu_pp, 96, 'filled', 'Color', color_cell_exvivo{1});
hold on;
plotHandles(:,2) = scatter(r2star_mean_illu_np, FF_mean_illu_np, 96, 'filled', 'Color', color_cell_avg16{1});
plotHandles(:,3) = scatter(r2star_mean_illu_nn, FF_mean_illu_nn, 96, 'filled', 'Color', color_cell_invivo{1});


set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_gray{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_gray{5}, 'MarkerFaceColor', color_cell_gray{2});

set(plotHandles(:,2), 'Marker', 's', 'Color', color_cell_exvivo{4});
set(plotHandles(:,2), 'Marker', 's',  'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_cell_exvivo{5}, 'MarkerFaceColor' , color_cell_exvivo{2});

set(plotHandles(:,3), 'Marker', '.', 'Color', color_cell_gray{4});
set(plotHandles(:,3), 'Marker', '^', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_gray{5}, 'MarkerFaceColor', color_cell_gray{2});

ylim([0 20]); xlim([0 100]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_np, FF_mean_illu_np, 1);
[r, p] = corrcoef(r2star_mean_illu_np, FF_mean_illu_np); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_np, FF_mean_illu_np, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_np);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_np) - FF_mean_illu_np).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 36:1:100;
y = x .* stats.slope + stats.intercept;
%plot(x, y, 'k', 'LineWidth', 2);
% axis off;

%% Purple
color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(r2star_mean_illu_pp, FF_mean_illu_pp, 96, 'filled', 'Color', color_cell_exvivo{1});
hold on;
plotHandles(:,2) = scatter(r2star_mean_illu_np, FF_mean_illu_np, 96, 'filled', 'Color', color_cell_avg16{1});
plotHandles(:,3) = scatter(r2star_mean_illu_nn, FF_mean_illu_nn, 96, 'filled', 'Color', color_cell_invivo{1});


set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_gray{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_gray{5}, 'MarkerFaceColor', color_cell_gray{2});

set(plotHandles(:,2), 'Marker', 's', 'Color', color_cell_gray{4});
set(plotHandles(:,2), 'Marker', 's',  'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_cell_gray{5}, 'MarkerFaceColor' , color_cell_gray{2});

set(plotHandles(:,3), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,3), 'Marker', '^', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

ylim([0 20]); xlim([0 100]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_nn, FF_mean_illu_nn, 1);
[r, p] = corrcoef(r2star_mean_illu_nn, FF_mean_illu_nn); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_nn, FF_mean_illu_nn, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_nn);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_nn) - FF_mean_illu_nn).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

%x = 0:0.2:4.3;
x = 0:1:34;
y = x .* stats.slope + stats.intercept;
%plot(x, y, 'k', 'LineWidth', 2);
%axis off;

%% Red T1 entropy vs FF
plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(FF_mean_illu_pp_v2, T1_entropy_illu_pp, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
xlim([0 20]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(FF_mean_illu_pp_v2, T1_entropy_illu_pp, 1);
[r, p] = corrcoef(FF_mean_illu_pp_v2, T1_entropy_illu_pp); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(FF_mean_illu_pp_v2, T1_entropy_illu_pp, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_pp_v2);
stats.SSE = sum((polyval(stats.polyCoefs, FF_mean_illu_pp_v2) - T1_entropy_illu_pp).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:20;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
axis square;
%% Red T1 entropy vs R2star

plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(r2star_mean_illu_pp_v2, T1_entropy_illu_pp, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
xlim([0 80]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_pp_v2, T1_entropy_illu_pp, 1);
[r, p] = corrcoef(r2star_mean_illu_pp_v2, T1_entropy_illu_pp); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_pp_v2, T1_entropy_illu_pp, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(r2star_mean_illu_pp_v2);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_pp_v2) - T1_entropy_illu_pp).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:80;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
axis square;

%% Blue T1 entropy vs FF


plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(FF_mean_illu_np_v2, T1_entropy_illu_np, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_exvivo{4});
set(plotHandles(:,1), 'Marker', 's', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_exvivo{5}, 'MarkerFaceColor', color_cell_exvivo{2});
xlim([0 20]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(FF_mean_illu_np_v2, T1_entropy_illu_np, 1);
[r, p] = corrcoef(FF_mean_illu_np_v2, T1_entropy_illu_np); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(FF_mean_illu_np_v2, T1_entropy_illu_np, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_np_v2);
stats.SSE = sum((polyval(stats.polyCoefs, FF_mean_illu_np_v2) - T1_entropy_illu_np).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:20;
y = x .* stats.slope + stats.intercept;
%plot(x,y, 'k', 'LineWidth', 2);
axis square;

%% Blue T1 entropy vs R2star
plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(r2star_mean_illu_np_v2, T1_entropy_illu_np, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_exvivo{4});
set(plotHandles(:,1), 'Marker', 's', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_exvivo{5}, 'MarkerFaceColor', color_cell_exvivo{2});
xlim([0 80]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_np_v2, T1_entropy_illu_np, 1);
[r, p] = corrcoef(r2star_mean_illu_np_v2, T1_entropy_illu_np); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_np_v2, T1_entropy_illu_np, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(r2star_mean_illu_np_v2);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_np_v2) - T1_entropy_illu_np).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:80;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
axis square;

%% Purple T1 entropy vs FF
plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(FF_mean_illu_nn_v2, T1_entropy_illu_nn, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,1), 'Marker', '^', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});
xlim([0 20]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(FF_mean_illu_nn_v2, T1_entropy_illu_nn, 1);
[r, p] = corrcoef(FF_mean_illu_nn_v2, T1_entropy_illu_nn); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(FF_mean_illu_nn_v2, T1_entropy_illu_nn, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_nn_v2);
stats.SSE = sum((polyval(stats.polyCoefs, FF_mean_illu_nn_v2) - T1_entropy_illu_nn).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:20;
y = x .* stats.slope + stats.intercept;
%plot(x,y, 'k', 'LineWidth', 2);
axis square;

%% Purple T1 entropy vs R2star
plotHandles = zeros(4,1);
figure();
plotHandles(:,1) = scatter(r2star_mean_illu_nn_v2, T1_entropy_illu_nn, 96);
set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,1), 'Marker', '^', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});
xlim([0 80]); ylim([1 5]);
hold on;

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu_nn_v2, T1_entropy_illu_nn, 1);
[r, p] = corrcoef(r2star_mean_illu_nn_v2, T1_entropy_illu_nn); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu_nn_v2, T1_entropy_illu_nn, 'type', 'Spearman');
stats.corrP = p(1,2);
stats.N = length(r2star_mean_illu_nn_v2);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu_nn_v2) - T1_entropy_illu_nn).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:80;
y = x .* stats.slope + stats.intercept;
%plot(x,y, 'k', 'LineWidth', 2);
axis square;

%% Examples Merry Slice 2
Merry_ff = [2.4998, 3.6279, 14.78];
Merry_r2star = [58.7306, 64.2459, 61.44];

color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(Merry_r2star, Merry_ff,  300, 'filled', 'Color', color_cell_exvivo{1});

set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

ylim([0 16]); xlim([0 70]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(Merry_r2star, Merry_ff, 1);
[r, p] = corrcoef(Merry_r2star, Merry_ff); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(Merry_r2star, Merry_ff);
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_pp);
stats.SSE = sum((polyval(stats.polyCoefs, Merry_r2star) - Merry_ff).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

%% Examples 18D15 Slice 2
ff_18D15 = [1.2185, 1.2974, 3.62];
r2star_18D15 = [45.2272, 61.9548, 47.87];

color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);

plotHandles(:,1) = scatter(r2star_18D15, ff_18D15,  300, 'filled', 'Color', color_cell_exvivo{1});

set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_exvivo{4});
set(plotHandles(:,1), 'Marker', 's', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_exvivo{5}, 'MarkerFaceColor', color_cell_exvivo{2});

ylim([0 16]); xlim([0 70]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_18D15, ff_18D15, 1);
[r, p] = corrcoef(r2star_18D15, ff_18D15); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_18D15, ff_18D15);
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_pp);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_18D15) - ff_18D15).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

%% Examples 18D16 Slice 2
ff_18D16 = [0.9345, 1.1306, 1.30];
r2star_18D16 = [18.1157, 29.6957, 30.2999];

color_cell_gray = {[245, 245, 245]/255, [220, 220, 220]/255, [211, 211, 211]/255, [192, 192, 192]/255, [169, 169 ,169]/255};
plotHandles = zeros(4,3);
%figure();
figure('Position', [0 100 600 550]);

plotHandles(:,1) = scatter(r2star_18D16, ff_18D16,  300, 'filled', 'Color', color_cell_invivo{1});

set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,1), 'Marker', '^', 'LineWidth', 2, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

ylim([0 16]); xlim([0 70]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_18D16, ff_18D16, 1);
[r, p] = corrcoef(r2star_18D16, ff_18D16); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_18D16, ff_18D16);
stats.corrP = p(1,2);
stats.N = length(FF_mean_illu_pp);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_18D16) - ff_18D16).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

%% Patient Analysis

% Remove 2 pixels
% T1_entropy = [3.226650665	3.457850828	2.718079221	2.499056589	3.371479193	3.195670399	3.064785548	2.804755242	2.975664237	3.392783749	3.377906439	3.008738456	2.98838875	2.844429901	2.832724493	2.697805734	2.822762112	3.029865506	2.81136788	2.954446942	2.725400327	2.409891914	2.474903441	2.778625884	2.656873042	2.252382635	2.70835284	3.030668964	2.639595553	2.359077838	3.071203999	2.84446803	2.436041881	2.763742186	1.992180025	2.203703282	1.993940438	1.785305929	1.986330983	2.02325357	1.998485279	2.457795368	2.290983872	2.198052043	2.581775734	2.678328269	2.115803637	1.332820405	2.516182382	2.625713045	2.010239527	1.914737636	3.197997566	3.505325142	3.074497691	3.089220972	2.908192061	2.84738103	2.937687855	2.596840402	3.000211946	3.169069808	2.994413901	2.524967107	2.998447499	2.565253142	2.20212135	1.835902812	2.329831286	2.427465917	2.198021807	1.87599672	3.175221734	3.322749077	2.377708927	3.182048804	2.95398955	2.545954336	2.951999651	3.290760029	3.080862045	2.380860235	3.293971748	3.007071329	2.189238256	2.485228136	3.241956654	2.850966975	2.509007051	1.557431962	3.330939307	3.090225788	3.0843122	3.30905914	3.409117229	2.773133387	1.959298911	0	1.554585169];

% Remove 1 pixel
% T1_entropy = [3.292015467	3.521195067	2.79463277	2.499056589	3.445044499	3.621144808	3.297626439	2.892417091	3.066471628	3.462203302	3.592682664	3.177138664	3.0530908	2.897942434	2.896518766	2.697805734	2.822762112	3.085583319	2.92440905	3.019313842	2.947289192	2.487492159	2.979806954	3.133844038	2.753393644	2.505249493	2.474903441	2.778625884	2.812441399	2.815071955	3.284571717	3.049826019	2.737397331	3.024822002	1.992180025	2.203703282	1.993940438	1.785305929	1.986330983	2.111291841	1.998485279	2.572645343	2.290983872	2.198052043	2.758429824	2.932287882	2.115803637	1.332820405	2.516182382	2.691009153	2.010239527	2.105397711	3.353331226	3.701725673	3.228027809	3.203343525	3.465365569	3.524141992	3.210316365	3.018648943	3.405099327	3.321758834	3.099583357	2.59655574	2.998447499	2.652450464	2.373915237	1.835902812	2.329831286	2.427465917	2.436638708	1.87599672	3.175221734	3.41357256	2.694810802	3.267397732	3.105601702	2.641704474	3.344269698	3.454493031	3.199919752	2.380860235	3.376717816	3.097800573	2.425098369	2.672905595	3.373857357	2.953169881	2.980240673	2.127065631	3.449855635	3.235358758	3.349767401	3.374464205	3.607389089	2.928387789	1.959298911	0.970950594	1.924174352];

% unweighted, 40 bins
% T1_entropy = [3.537273219	4.032394521	3.326644545	3.008024795	3.927864033	4.258540624	3.857469858	3.479601779	3.799957565	4.233525404	4.002220303	3.694694354	3.571369543	3.598695476	3.482047247	3.367627438	3.436607051	3.43726271	3.513891595	3.47255582	3.540584069	3.085670113	3.288293569	3.593676212	3.259061911	2.704099781	2.617536533	3.276420835	3.067476457	2.652345561	3.928134475	3.478835254	3.43127238	3.768273847	2.4758585	2.74216951	2.303627993	2.232848992	2.383988609	2.682792823	2.516696304	3.079797834	2.906915359	2.826764657	3.365261528	3.604876133	2.936911074	2.214199695	3.118697829	3.342514592	2.697405766	2.556359143	3.502596841	4.049195895	3.654991497	3.415678079	3.694891829	4.018089708	3.975646384	3.30708252	4.047772969	3.940762991	3.723397038	3.237951205	3.56244661	3.198341371	2.739881706	2.271066057	2.561964212	2.893388358	2.879451104	2.765574736	3.190119875	3.52847719	3.197407646	3.392662215	3.714523193	3.109386672	3.564926018	4.125803243	3.501826416	2.928782852	4.001909518	3.65633558	3.374548236	2.96147252	3.995461119	3.552783151	3.358922883	3.269576746	3.66462172	3.290420494	3.343119911	3.55559154	3.427307576	3.022182094	2.307820398	2.05881389	3.090032777];

% unweighted remove 1 pixels
% T1_entropy = [3.537273219	3.782380131	3.286405173	2.888662706	3.84755692	3.864413376	3.777840376	3.432148197	3.366493661	3.944728716	3.725613275	3.458720552	3.394581786	3.484402087	3.447749151	3.226404779	3.352014857	3.310557034	3.394060013	3.436646427	3.463385702	3.044935395	3.137126814	3.294681477	2.934717607	2.532314007	2.36559623	3.029408653	2.772208352	2.41268789	3.692190757	3.421283068	3.356971263	3.56895206	2.288466328	2.529550926	2.220323638	2.15218353	2.383988609	2.588147006	2.446481784	3.014843837	2.781837408	2.785245314	3.273114274	3.604876133	2.796907293	2.086641775	2.962271231	3.308123598	2.697405766	2.236779563	3.3193688	3.898256069	3.570623578	3.131673258	3.239668316	3.617076226	3.840960127	3.074820572	3.896728115	3.73576803	3.494882387	3.08415993	3.486437276	3.103631795	2.64618737	2.271066057	2.354435538	2.83859057	2.751154309	2.560352492	3.065543713	3.360131473	3.030869578	3.340515087	3.388982986	2.898502483	3.324424043	3.818738355	3.140864629	2.627413428	3.824450462	3.554838897	2.503997527	2.727827978	3.89027626	3.327108655	3.103347507	2.935649892	3.525242792	3.247439407	3.303235567	3.479055348	3.427307576	2.975447368	1.94400975	0.918295834	2.235926351];

% Unweighted, 20 bins

T1_entropy = [2.28205291	2.50214416	2.69884046	3.56599615	2.10720922	2.51740581	2.89180359	2.68755764	2.07183696	2.43498128	2.66059441	3.0554517	3.25769513	3.16850826	3.14854754	3.05574655	2.76399946	3.06144447	3.0766157	3.21425814	2.72657241	2.68947245	3.4038486	3.49321462	2.32291632	2.44754466	3.0929984	3.51550082	3.17677997	2.84691924	2.81413748	2.88697022	2.8535377	2.3216567	3.52521305	3.47485847	2.97466235	3.092948	2.9045812	3.06398577	2.58429224	3.18021263	3.08652287	3.11695612	2.00118504	2.83910566	2.36984604	2.55803844	2.29338114	2.81278227	3.11698756	3.1791488	2.58384232	2.33675755	3.13986922	3.32012316	2.53483908	3.73822905	2.48163653	3.06092667	3.02712856	3.34356331	3.2417901	3.08028296	3.21981438	2.50558227	3.28995237	2.82922865	3.1724382	2.4972527	2.56598725	2.81009045	2.86585417	2.45111059	2.9032683	3.02932544	3.69591776	3.64955774	3.04907748	2.72255453	3.18680658	3.13742148	2.68536749	2.49959198	2.75658587	2.81830747	3.13528262	3.05221568];
%T1_entropy_remote = [0.865856617	0.67694187	1.480290753	1.2586633	1.236386411	1.423794941	1.681046078	1.251629167	1.494918848	1.337778097	1.700438513	1.323489961	1.58468025	1.672492724	1.429152984	2.134711144	1.952819531	1.334599426	1.620538846	1.901872151	1.54536616	1.846393983	1.987773371	1.945661003	2.022887961	1.294995075	1.422475118	1.573114618	1.916445381	1.040852083	1.481863125	1.696842803	1.870018047	0.985228136	0.99631652	1.253297578	2.230337797	1.530639062	1.086312841	0.976020648	1.392147224	2.43147771	1.726474118	1.231724487	0.811847521	1.485475297	1.021119189	0.811278124	1.337306179	0.73206669	1.298794941	1.052981575	0.937301241	0.964086328	1.303169878	1.224394445	0.787004656	1.358584608	1.632105718	1.376074607	1.924820746	2.095646899	1.217016339	0.793992434	1.030426461	0.970950594	1.56909404	0.959686894	1.722971279	1.226817429	1.79650626	1.432928309	2.09790657	1.571091977	1.742966195	1.646780976	1.67673703	1.214540408	1.468590448	1.254999238	1.130329644	1.882045108	1.217038755	1.109227139	1.115644407	1.469182539	1.570950594	1.421554169	1.248609233	2.087210311	1.781839569	1.906305535	0.991264261	1.514723983	0.959686894	1.325784781	1.498395816	1.467457965	0.266764988	1.582683189];

FF_mean = [2.38395067	2.46726384	1.39470268	4.90014523	2.4	2.68491684	7.55650747	7.62215284	1.4228047	3.1842703	3.70914201	8.02984403	8.8982707	9.783539	7.19338363	3.37815015	2.58447272	2.64623091	2.8487631	4.47691648	3.2279974	3.01563118	2.53382371	3.01823052	1.75416354	2.73620219	3.66549854	4.18963012	4.99352532	3.53700525	1.65291897	2.30068687	2.83098455	2.94785235	13.6921508	13.8069145	6.69161525	5.60829414	1.19072113	3.1337654	3.107209	3.6576403	3.46242685	2.40637669	2.81210708	2.05996528	2.05227734	4.06120379	4.28592579	3.79310693	2.30706634	3.60990906	1.29560677	2.48523145	1.66743548	5.7046877	1.63289635	2.42853373	2.12635934	3.4687958	2.14086753	3.04401981	3.08940503	5.30876759	2.84262002	3.36181338	3.33375328	5.7530123	2.77157061	3.54290923	2.73971731	3.87800649	1.84170094	2.62871486	3.07179554	3.46324735	4.15715748	3.44925921	2.65575733	3.02652666	1.48482578	2.69588627	5.14265604	3.65985824	3.28511144	2.79396432	8.88611321	4.32083933];

r2star_mean = [29.7020392	42.100792	28.0536824	32.1219663	32.3716792	35.0245622	45.1212035	37.7591537	32.2862123	35.3344049	40.4318816	49.8894561	45.3424275	42.8906535	36.4294489	31.7247566	27.3249777	26.8516967	26.0167722	29.685374	32.2657114	32.2131726	25.6647431	27.6616562	32.7849176	31.6740012	35.6561527	30.6249763	39.1911374	46.8802944	18.4349467	26.2321835	28.87341	18.1768863	59.5624289	61.0137078	43.0116781	35.04743	12.6183619	31.1019292	31.3651919	33.9852718	34.6382513	35.9416697	23.4992489	27.4311097	24.6292201	27.8068147	20.8040457	28.1670478	35.3129534	29.1575121	25.0938706	30.0745682	30.3482152	36.0100905	30.1987607	29.11512	28.0487211	41.5535327	41.2484596	32.2687617	29.6705669	48.1701859	36.1075937	43.4147908	53.2384405	31.291723	32.0686488	32.4730604	31.1359025	29.1417872	27.1384095	34.8471951	42.2322058	38.7603209	49.3848695	40.0680296	35.6858058	32.1516574	20.8021	23.8998422	45.272934	31.5279943	29.9504749	34.2141938	42.1676878	50.54663];

find(abs(zscore(T1_entropy))>3)
find(abs(zscore(FF_mean))>3)     % 7, 8
find(abs(zscore(r2star_mean))>3) % 79, 80
% find(abs(zscore(LGE_entropy))>3) % 87, 98
% find(abs(zscore(T2_mean))>3) % 50

find(abs(zscore(T1_entropy))>2.5)  % 48
find(abs(zscore(FF_mean))>2.5)     % 7, 8, 9
find(abs(zscore(r2star_mean))>2.5) % 79, 80
% find(abs(zscore(LGE_entropy))>2.5) % 87, 98
% find(abs(zscore(T2_mean))>2.5) % 50

find(abs(zscore(T1_entropy))>2)  % 38, 48, 68
find(abs(zscore(FF_mean))>2)     % 6, 7, 8, 9, 10
find(abs(zscore(r2star_mean))>2) % 31, 79, 80
% find(abs(zscore(LGE_entropy))>2) % 33, 87, 98
% find(abs(zscore(T2_mean))>2) % 49, 50, 51, 78, 87

exclude_idx_ff = zeros(1, length(FF_mean)) > 0;
exclude_idx_r2star = zeros(1, length(r2star_mean)) > 0;
exclude_idx = zeros(1, length(FF_mean)) > 0;


% exclude_idx_ff = FF_mean < 4;
% exclude_idx_r2star = r2star_mean < 35; % T1_entropy = 2.5968
% 
% exclude_idx_ff = FF_mean >= 4;
% exclude_idx_r2star = r2star_mean >= 35; % T1_entropy = 2.1776

% exclude_idx_ff = FF_mean < 4;
% exclude_idx_r2star = r2star_mean >= 35;
% 
% exclude_idx_ff = FF_mean >= 4;
% exclude_idx_r2star = r2star_mean < 35; % T1_entropy = 2.1552

% exclude_idx_r2star = r2star_mean < 40;
exclude_idx_hemo = hemo_label;
exclude_idx = exclude_idx_ff | exclude_idx_r2star;

exclude_idx(77) = 1;
% exclude_idx(82) = 1;

T1_entropy(exclude_idx) = [];
r2star_mean(exclude_idx) = [];
FF_mean(exclude_idx) = [];

%% For publication (Patient)
color_cell_exvivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
color_cell_invivo = {[242,240,247]/255, [203,201,226]/255, [158,154,200]/255, [117,107,177]/255, [84,39,143]/255}; % Purple
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};

FF_mean_illu = zeros(1, length(FF_mean));
FF_mean_illu_label_pp = find(FF_mean >= 4 & r2star_mean >= 35);
FF_mean_illu_label_np = find(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_label_nn = find(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu_pp = FF_mean(FF_mean >= 4 & r2star_mean >= 35) + 0.6;
FF_mean_illu_np = FF_mean(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_nn = FF_mean(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu = zeros(1, length(r2star_mean));
r2star_mean_illu_label_pp = find(FF_mean >= 4 & r2star_mean >= 35);
r2star_mean_illu_label_np = find(FF_mean < 4 & r2star_mean >= 35);
r2star_mean_illu_label_nn = find(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu_pp = r2star_mean(FF_mean >= 4 & r2star_mean >= 35) + 3;
r2star_mean_illu_np = r2star_mean(FF_mean < 4 & r2star_mean >= 35) + 3;
r2star_mean_illu_nn = r2star_mean(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu(FF_mean_illu_label_pp) = FF_mean_illu_pp;
FF_mean_illu(FF_mean_illu_label_np) = FF_mean_illu_np;
FF_mean_illu(FF_mean_illu_label_nn) = FF_mean_illu_nn;

r2star_mean_illu(r2star_mean_illu_label_pp) = r2star_mean_illu_pp;
r2star_mean_illu(r2star_mean_illu_label_np) = r2star_mean_illu_np;
r2star_mean_illu(r2star_mean_illu_label_nn) = r2star_mean_illu_nn;

figure();
scatter(FF_mean_illu, r2star_mean_illu, 'filled');
xlim([0 20]); ylim([0 100]);
hold on;
xline(4.25);
yline(36);

T1_entropy_illu_pp = T1_entropy(FF_mean >= 4 & r2star_mean >= 35);
T1_entropy_illu_np = T1_entropy(FF_mean < 4 & r2star_mean >= 35);
T1_entropy_illu_nn = T1_entropy(FF_mean < 4 & r2star_mean < 35);

% T2_mean_illu_pp = T2_mean(FF_mean >= 4 & r2star_mean >= 35);
% T2_mean_illu_np = T2_mean(FF_mean < 4 & r2star_mean >= 35);
% T2_mean_illu_nn = T2_mean(FF_mean < 4 & r2star_mean < 35);

FF_mean_illu_pp_v2 = FF_mean(FF_mean >= 4 & r2star_mean >= 35);
FF_mean_illu_np_v2 = FF_mean(FF_mean < 4 & r2star_mean >= 35);
FF_mean_illu_nn_v2 = FF_mean(FF_mean < 4 & r2star_mean < 35);

r2star_mean_illu_pp_v2 = r2star_mean(FF_mean >= 4 & r2star_mean >= 35);
r2star_mean_illu_np_v2 = r2star_mean(FF_mean < 4 & r2star_mean >= 35);
r2star_mean_illu_nn_v2 = r2star_mean(FF_mean < 4 & r2star_mean < 35);

% Exclusion
FF_mean_illu_np_v2([]) = [];
T1_entropy_illu_np([]) = [];
r2star_mean_illu_np_v2([]) = [];

FF_mean_illu_np([]) = [];
r2star_mean_illu_np([]) = [];

FF_mean_illu_pp_v2([14]) = [];
r2star_mean_illu_pp_v2([14]) = [];
T1_entropy_illu_pp([14]) = [];
%% For publication
plotHandles = zeros(4,3);
% figure();
figure('Position', [0 100 600 550]);
plotHandles(:,1) = scatter(r2star_mean_illu_pp, FF_mean_illu_pp,  64, 'filled', 'Color', color_cell_exvivo{1});
hold on;
plotHandles(:,2) = scatter(r2star_mean_illu_np, FF_mean_illu_np, 64, 'filled', 'Color', color_cell_avg16{1});
plotHandles(:,3) = scatter(r2star_mean_illu_nn, FF_mean_illu_nn, 64, 'filled', 'Color', color_cell_invivo{1});

set(plotHandles(:,2), 'Marker', 's', 'Color', color_cell_exvivo{4});
set(plotHandles(:,2), 'Marker', 's',  'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_cell_exvivo{5}, 'MarkerFaceColor' , color_cell_exvivo{2});


set(plotHandles(:,1), 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'Marker', 'o', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

set(plotHandles(:,3), 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,3), 'Marker', '^', 'LineWidth', 1.5, ...
    'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

ylim([0 20]); xlim([0 100]);
hold on;
yline(4.25);
xline(36);

stats = struct;
[stats.polyCoefs, stats.polyFitStruct] = polyfit(r2star_mean_illu, FF_mean_illu, 1);
[r, p] = corrcoef(r2star_mean_illu, FF_mean_illu); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(r2star_mean_illu, FF_mean_illu);
stats.corrP = p(1,2);
stats.N = length(r2star_mean_illu);
stats.SSE = sum((polyval(stats.polyCoefs, r2star_mean_illu) - FF_mean_illu).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);

x = 0:1:100;
y = x .* stats.slope + stats.intercept;
plot(x,y, 'k', 'LineWidth', 2);
%axis off;
%% Please go to line 414 %Red for more information
