%% LoopProgressReport
%==========================================================================
% This function prints the loop progress. Typical usage:
% LoopProgressReport('Staring loop')
% for i = 1:100
%   %Do something
%   LoopProgressReport(i,100)
% end
% It uses the tic/toc function in Matlab, this means that if tic is called
% inside the loop, the function is not going to work correctly. The time
% left estimate is based on the asumption that each itterations is equally
% long.
% Created by Anders, Jul 4. 2016
%==========================================================================

function   LoopProgressReport(itterNo,totalItters)
warning off
%initialization without custom string:
if nargin==0,		% initialisation
    [ST, I] = dbstack('-completenames', 1);
    fprintf([ST(1).name ': ' 'Loop execution: 0%%  ']);
    tic();
    Percent = 0;
end

%initialization with custom string:
if nargin==1,
    titleString = itterNo;
    %Check if input is string
    if ~ischar(titleString),
        error('Single input parameter must be a string');
    end
    
    %Add : at the end (but ensure that we dont get two ::)
    titleString = strrep(titleString,char(10),'');
    titleString = [strrep(titleString, ': ','') ': '];
    [ST, I] = dbstack('-completenames', 1);
    try
        fprintf([ST(1).name ': ' titleString '          ' ]);
    catch
        fprintf(['Script: ' titleString '          ' ]);
    end
    timeLeftValue =round(0/60);
    timeLeftString = ['                    ' num2str(timeLeftValue)];
    timeLeftString = timeLeftString(end-3:end);
    timeLeftString = ['           '  '      '];
    fprintf(timeLeftString)
    
    tic();
end

%For every itteration
if nargin==2,
    
    % use with FOR loop
    if isempty(getCurrentTask());
        
        Percent = floor(100*itterNo/totalItters);
        PercentOld = floor(100*(itterNo-1)/totalItters);
        if PercentOld~=Percent || Percent == 0 && itterNo > 1
            %Time to now:
            
            t = toc;
            
            %Time pr itter:
            t_ = toc/itterNo;
            
            %itters left?
            left = totalItters-itterNo;
            
            %Time left
            t = left*t_;
            
            %Hours
            number = floor(t/60/60);
            suffix = 'h  ';
            
            %minutes
            if number == 0
                number = floor(t/60);
                suffix = 'min';
                %seconds
                if number == 0
                    number = round(t);
                    suffix = 's  ';
                else
                    number = floor(t/60);
                end
                
            
            else
                number = round(t/60/60);
            end
            
            
            
            
            timeLeftString = ['                    ' num2str(number)];
            timeLeftString = timeLeftString(end-3:end);
            timeLeftString = [' Time left:'  timeLeftString suffix];
            backspace = '\b';
            
            percentString = ['         ' num2str(Percent) '%%'];
            percentString = percentString(end-4:end);
            
            fprintf([repmat(backspace, 1,length(timeLeftString)+length(percentString)-1) percentString ]);
            
            fprintf(timeLeftString)
            
            if itterNo==totalItters,	% safe to say that this is the last call to PrintLoopPC...
                %fprintf(char(10));
                
                t = toc;
                
                %Hours
                number = round(t/60/60);
                suffix = 'h';
                
                %minutes
                if number == 0
                    number = round(t/60);
                    suffix = 'min';
                end
                
                %seconds
                if number == 0
                    number = round(t);
                    suffix = 's';
                end
                
                
                fprintf([ repmat('\b',1,22) '  ' 'Finished after ' num2str(number) suffix  char(10)])
            end
            
            
        end
    else %If parfor loop do nothing..
        
    end
end
warning on
end
