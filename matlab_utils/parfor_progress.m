function percent = parfor_progress(pfad, N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(pfad, N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS(pfad) updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(pfad, 0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      pfad = 'c:\temp';
%      parfor_progress(pfad, N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress(pfad);
%      end
%      parfor_progress(pfad, 0);
%
%   See also PARFOR.


if nargin < 2 && ischar(pfad)
    narginchk(1, 1);
elseif nargin < 3 && ischar(pfad) && (isnumeric(N) || islogical(N))
    narginchk(2, 2);
else
    error('Falsche Übergabeparameter!');
end

replace = false;
if nargin < 2
    N = -1;
elseif islogical(N)
    replace = N;
    N = -1;
end

percent = 0;
w = 73; % Width of progress bar

if N > 0
    if ~exist(fullfile(pfad), 'dir')
        mkdir(pfad);
    end
    f = fopen(fullfile(pfad, 'parfor_progress.txt'), 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete(fullfile(pfad, 'parfor_progress.txt'));
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
elseif replace
    anzahl = 0;
    while exist(fullfile(pfad, 'parfor_progress.txt'), 'file') == 0
        if anzahl > 9
            fprintf(2, 'parfor_progress.txt not found.');
            return
        end
        anzahl = anzahl + 1;
        pause(0.3);
    end
    
    f = fopen(fullfile(pfad, 'parfor_progress.txt'), 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([ perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
else
    if ~exist(fullfile(pfad, 'parfor_progress.txt'), 'file')
        fprintf('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
        return
    end
    
    f = fopen(fullfile(pfad, 'parfor_progress.txt'), 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen(fullfile(pfad, 'parfor_progress.txt'), 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end
