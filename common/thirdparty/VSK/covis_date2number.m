function datevalue=covis_date2number(month,day,hour)
%
% Transform the input date (month/day/hour) into a numerical value in terms
% of the number of hours each date experienced from the beginning of a year
%
% ---------------------------
% written by Li Liu in 01/11/2013 
% l.liu6819@gmail.com
%

%compute the hours of each month
switch month
    case 1,
        m=0;
    case 2,
        m=24*31;
    case 3,
        m=24*59;
    case 4,
        m=24*90;
    case 5,
        m=24*120;
    case 6,
        m=24*151;
    case 7,
        m=24*181;
    case 8,
        m=24*212;
    case 9,
        m=24*243;
    case 10,
        m=24*273;
    case 11,
        m=24*304;
    case 12,
        m=24*334;
end

%compute the hours of every day
d=24*day;

datevalue=m+d+hour;
