function Prediction_chart(Data_strain)
ORF_name = Data_strain{1}; % ORF_name_Column A
Gene_name = Data_strain{2}; % Gene_name_Column B
Seq_lenght = Data_strain{4}; % Seq_lenght_Column D 
signal_p = Data_strain{14}; % signal_peptide_lengh TargetP(SP) column N
Probability_signal_p = Data_strain{15}; % Probability of TargetP (SP) column O
Validation_for_Nterminal = Data_strain{17}; % Experimental validation of N’ “IN” orientation  column Q
Validation_for_MTS = Data_strain{20}; % Experimental validation of MTS column T
Probability_target_MTS = Data_strain{21}; % Probability TargetP (SP) column U
lenght_target_MTS = Data_strain{22}; % Lenght TargetP (SP) column V
Probability_MitoFates = Data_strain{23}; % Probability MTS column W
Lenght_MitoFates = Data_strain{24}; % Lenght MTS column X
Description = Data_strain{25}; % Description column Y

Matrix_predictions = {}; % Rows=type of prediction, Coloumn=Location of the characters, when 1=M 2=S 3=i 4=o.
for k=1:7    
    curent_prediction = Data_strain{6+k};
    Matrix_predictions{k,1} = strfind(curent_prediction,'M');
    Matrix_predictions{k,2} = strfind(curent_prediction,'S');
    Matrix_predictions{k,3} = strfind(curent_prediction,'i');
    Matrix_predictions{k,4} = strfind(curent_prediction,'o');
end

f = figure('Position',[337 20 1248 1096]);
name_of_predictions = {'TMHMM','TOPCONS','OCTOPUS','Philius','PolyPhobius','SCAMPI','SPOCTOPUS'};
cordination_for_y_axis = {[0.918012412668042 0.882921384309793 0.0599190267956692 0.0260617755727418],...
    [0.914345138233556 0.799005444918713 0.0712550588464931 0.0260617755727418],...
	[0.910140350717996 0.718887629710924 0.0712550588464931 0.0260617755727418],...
	[0.912597314743863 0.637930616877305 0.0534412941951984 0.0260617755727418],...
	[0.904998429464832 0.552313163192832 0.0834008074723758 0.0260617755727418],...
	[0.911487279750265 0.471863400151639 0.0599190267956692 0.0260617755727418],...
	[0.909867846551883 0.384990813279052 0.0858299571975524 0.0260617755727418],};
%% 6 Predictions programs                    
for l = 1:7
    plot_seq = subplot(10,1,l);
    M = Matrix_predictions{l,1};
    S = Matrix_predictions{l,2};
    i = Matrix_predictions{l,3};
    o = Matrix_predictions{l,4};
    scatter(M,ones(1,length(M))*1,200,[0.560784339904785 0.819607853889465 0.603921592235565],'s','filled')
    hold on
    scatter(S,ones(1,length(S))*1,200,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
    scatter(i,ones(1,length(i))*1,200,[0.929411768913269 0.686274528503418 0.686274528503418],'s','filled')
    scatter(o,ones(1,length(o))*1,200,[0.603921592235565 0.603921592235565 0.792156875133514],'s','filled')
    
    Change_points = [];
    for in_sqe = 2:Seq_lenght
        if Data_strain{6+l}(in_sqe) ~= Data_strain{6+l}(in_sqe-1)
            Change_points = [Change_points,in_sqe];
        end
    end
    Change_points = [Change_points];
    
    annotation('textbox',cordination_for_y_axis{l},'String',name_of_predictions{l},'FitBoxToText','on','EdgeColor','none','FontWeight','bold');
    box on
    range_seq = [-Seq_lenght/100 Seq_lenght+Seq_lenght/100];
    xlim(range_seq);
    ylim([0,1.5])
    if l == 7
        set(gca,'xgrid','off');
    else
        set(gca,'xtick',[]);
    end
    
    set(gca,'ytick',[]);
    
    if not(isempty(Change_points))
        stem(Change_points-1,ones(1,length(Change_points))*1,'Color','r','Marker','o','LineStyle','-.',...
            'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
        state = 0.2;
        for in_sqe = 1:length(Change_points)
            if in_sqe == 1
                text_level = 0.2;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','bold','Color','b');
                continue
            end
            
            if Change_points(in_sqe)-Change_points(in_sqe-1) < 25 && state ~= 0.4
                text_level = 0.4;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','bold','Color','b');
                state = 0.4;
            else 
                text_level = 0.2;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','bold','Color','b');
                state = 0.2;
            end     
        end
    end
end

Original_tick = plot_seq.XTick;
if Seq_lenght <= Original_tick(end)        
    plot_seq.XTick(end) = Seq_lenght;
else
    plot_seq.XTick = [plot_seq.XTick,Seq_lenght];
end
    if (plot_seq.XTick(end)-plot_seq.XTick(end-1)) < 20
    plot_seq.XTick(end-1) = [];
end

% signal_p
   cur_plot = subplot(10,1,8);
   cur_pos = cur_plot.Position;
   set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.4,cur_pos(4)]);
if signal_p ~= 0
   scatter(1:signal_p,ones(1,signal_p)*1,200,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
   hold on
   stem(signal_p,1,'Color','r','Marker','o','LineStyle','-.',...
   'LineWidth',0.8,'MarkerSize',8,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
   text(signal_p,0.2,num2str(signal_p),'FontWeight','bold','Color','b');
   % scatter(signal_p,ones(1,length(signal_p))*1,200,'y','s','filled')
end

annotation('textbox',[0.516859975305978 0.269305019305019 0.0872745030123045 0.0658057984536131],...
    'String',{'TargetP (SP)',['probability: ',num2str(Probability_signal_p)]},...
    'EdgeColor','none','FontWeight','bold');
set(gca,'ytick',[])
xlim([-3,100]);
ylim([0,1.5]);
set(gca,'xtick',[]);
box on
% probabilities target, MTS
subplot(10,1,9)
cur_plot=subplot(10,1,9);
cur_pos=cur_plot.Position;
set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.4,cur_pos(4)])

if lenght_target_MTS ~= 0
    scatter(1:lenght_target_MTS,ones(1,lenght_target_MTS)*1,200,[0.929411764705882 0.690196078431373 0.129411764705882],'s','filled')
    hold on
    stem(lenght_target_MTS,1,'Color','r','Marker','o','LineStyle','-.',...
    'LineWidth',0.8,'MarkerSize',8,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
    text(lenght_target_MTS,0.2,num2str(lenght_target_MTS),'FontWeight','bold','Color','b');
end

annotation('textbox',...
    [0.516859975305978 0.180501930501931 0.0894970697142401 0.0748231382318783],'String',{'TargetP (MTS)','probability:',num2str(Probability_target_MTS)},...
    'EdgeColor','none','FontWeight','bold');
set(gca,'ytick',[])
xlim([-3,100]);
ylim([0,1.5]);
set(gca,'xtick',[]);
box on
subplot(10,1,10)
cur_plot=subplot(10,1,10);
cur_pos=cur_plot.Position;
set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.4,cur_pos(4)])

if Probability_MitoFates ~= 0
    scatter(1:Lenght_MitoFates,ones(1,Lenght_MitoFates)*1,200,[0.929411764705882 0.690196078431373 0.129411764705882],'s','filled')
end

set(gca,'ytick',[])
xlim(range_seq);
hold on

if Probability_MitoFates ~= 0
    stem(Lenght_MitoFates,1,'Color','r','Marker','o','LineStyle','-.',...
        'LineWidth',0.8,'MarkerSize',8,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
    text(Lenght_MitoFates,0.2,num2str(Lenght_MitoFates),'FontWeight','bold','Color','b');
end

box on
xlim([-3,100]);
ylim([0,1.5]);  
annotation('textbox',[0.516859975305978 0.0777262756697201 0.0495369751021926 0.0914498714652955],...
    'String',{'MitoFates',['probability: ',num2str(Probability_MitoFates)]},...
    'EdgeColor','none','FontWeight','bold');
% legend
hold on
annotation('textbox',...
    [0.631873058277303 0.271218151680877 0.281486453408256 0.0504347351391072],'Color',[0.603921592235565 0.603921592235565 0.792156875133514],...
    'String','Out (Luminal or Extra-Cellular)',...
    'FitBoxToText','off','EdgeColor','k',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);

annotation('textbox',...
    [0.631873058277303 0.260555726494028 0.306654820350932 0.0397077349108242],...
    'Color',[0.560784339904785 0.819607853889465 0.603921592235565],...
    'String','Trans Membrane Domain',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);

annotation('textbox',[0.631873058277303 0.244526108921995 0.29460573683777 0.0321686931147077],...
    'Color',[0.729411780834198 0.831372559070587 0.95686274766922],...
    'String','Cleavable Signal peptide (SP)',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);

annotation('textbox',...
    [0.631873058277303 0.316170559101151 0.294544563704537 0.0278971027582894],...
    'Color',[0.929411768913269 0.686274528503418 0.686274528503418],...
    'String','In (Cytosol)',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);

annotation('textbox',...
    [0.631873058277303 0.214973052376651 0.275567170580216 0.0382111443062618],...
    'Color',[0.929411764705882 0.690196078431373 0.129411764705882],...
    'String','Mitochondrial Targeting Sequence (MTS) ',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
%% the rest of the data
% title
hold on
annotation('textbox',...
    [0.131668119324453 0.950772200772203 0.124753047290049 0.049637065637071],...
    'String',{[ORF_name,' ',Gene_name]},...
    'EdgeColor','none','FontWeight','bold','FontSize',17);      


annotation('textbox',...
    [0.251653423859989 0.944980694980695 0.644968981471828 0.0524285714285755],...
    'String',{Description},...
    'EdgeColor','none');

% expermintal validation Results: 4-NO, 3-vogtle 2009, 2-Venna 2013, 1-Yes both
% Experimental validation for n'
if Validation_for_MTS == 1
    annotation('textbox',[0.631873058277303 0.166023166023167 0.283120485249518 0.0558687258687266],...
    'String',{'Experimental validation of MTS: Yes, Venne2013, Vogtle 2009'},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
elseif Validation_for_MTS == 2
        annotation('textbox',[0.631873058277303 0.166023166023167 0.283120485249518 0.0558687258687266],...
    'String',{'Experimental validation of MTS: Yes, Vogtle 2009'},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
elseif Validation_for_MTS == 3
        annotation('textbox',[0.631873058277303 0.166023166023167 0.283120485249518 0.0558687258687266],...
    'String',{'Experimental validation of MTS: Yes,  Venne2013'},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
elseif Validation_for_MTS == 4
        annotation('textbox',[0.631873058277303 0.166023166023167 0.283120485249518 0.0558687258687266],...
    'String',{'Experimental validation of MTS: No'},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
end
    
if Validation_for_Nterminal == 3
    annotation('textbox',...
    [0.631873058277303 0.138482625482627 0.341777777777778 0.0300000000000002],'String',{"Experimental validation of N’ “IN” orientation: Not available"},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
elseif Validation_for_Nterminal == 2
    annotation('textbox',...
    [0.631873058277303 0.138482625482627 0.341777777777778 0.0300000000000002],'String',{"Experimental validation of N’ “IN” orientation: Yes"},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
elseif Validation_for_Nterminal == 1
        annotation('textbox',...
    [0.631873058277303 0.138482625482627 0.341777777777778 0.0300000000000002],'String',{"Experimental validation of N’ “IN” orientation: No"},...
    'EdgeColor','none','FontWeight','bold','FontSize',13);
end