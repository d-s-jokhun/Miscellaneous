%%%written by D.S.Jokhun on 27.07.2017 %%%
%%% finds out within what range can a graph be considered linear


dataX=timescale(10:end,:);
% dataY=SD(10:end,:);

SD_MA=[];
for count=1:size(SD,2)
SD_MA(:,count)=smooth(SD(:,count),5,'moving');
end
dataY=SD_MA(10:end,:);





%% just need to make sure 'x' and 'y' in the section below contains the curve to be analysed at every loop

inflection_criterion1 = 1e-4;  % sse value beyond which the linear fit will not be considered good, meaning that after this point, it's not a line but a curve.
inflection_criterion2 = 5e-4;

result_sse=[];
result_inflection_point=[];  %in terms of x
for count_cond=1:size(dataY,2)
    check=0;
    for count_time=2:size(dataX,1)
        x=[];
        y=[];
        x=dataX(1:count_time,1);
        y=dataY(1:count_time,count_cond);
        [fitobject,gof] = fit(x,y,'poly1');
        result_sse(count_time,count_cond)=gof.sse;
        
        if check==0
            if gof.sse >= inflection_criterion1
                check=1;
                result_inflection_point(1,count_cond)=x(count_time,1);
            end
        end
        
        if check==1
            if gof.sse >= inflection_criterion2
                check=2;
                result_inflection_point(2,count_cond)=x(count_time,1);
            end
        end
        
        
    end
end











