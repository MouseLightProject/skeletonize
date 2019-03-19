function date=covis_number2date(datevalue)
%
% Transform the input datevalue into date(month/day).
%
% ---------------------------
% written by Li Liu in 01/11/2013 
% l.liu6819@gmail.com
%

hour=mod(datevalue,24);

%compute the corresponding integer number of days of the input datevalue
totalday=(datevalue-hour)/24;

%compute month of the datevalue according to number of days each month has
if totalday<31,
    day=totalday;
    Month=strcat(num2str(0),num2str(1));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==31
    Month=strcat(num2str(0),num2str(1));
    Day=num2str(31);
end
    
if totalday>31 && totalday<59,
    day=totalday-31;
    Month=strcat(num2str(0),num2str(2));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==59
    Month=strcat(num2str(0),num2str(2));
    Day=num2str(28);
end

if totalday>59 && totalday<90,
    day=totalday-59;
    Month=strcat(num2str(0),num2str(3));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==90
    Month=strcat(num2str(0),num2str(3));
    Day=num2str(31);
end

if totalday>90 && totalday<120,
    day=totalday-90;
    Month=strcat(num2str(0),num2str(4));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==120
    Month=strcat(num2str(0),num2str(4));
    Day=num2str(30);
end

if totalday>120 && totalday<151,
    day=totalday-120;
    Month=strcat(num2str(0),num2str(5));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==151
    Month=strcat(num2str(0),num2str(5));
    Day=num2str(31);
end

if totalday>151 && totalday<181,
    day=totalday-151;
    Month=strcat(num2str(0),num2str(6));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==181
    Month=strcat(num2str(0),num2str(6));
    Day=num2str(30);
end

if totalday>181 && totalday<212,
    day=totalday-181;
    Month=strcat(num2str(0),num2str(7));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==212
    Month=strcat(num2str(0),num2str(7));
    Day=num2str(31);
end

if totalday>212 && totalday<243,
    day=totalday-212;
    Month=strcat(num2str(0),num2str(8));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==243
    Month=strcat(num2str(0),num2str(8));
    Day=num2str(31);
end

if totalday>243 && totalday<273,
    day=totalday-243;
    Month=strcat(num2str(0),num2str(9));
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==273
    Month=strcat(num2str(0),num2str(9));
    Day=num2str(30);
end

if totalday>273 && totalday<304,
    day=totalday-273;
    Month=num2str(10);
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==304
    Month=num2str(10);
    Day=num2str(31);
end

if totalday>304 && totalday<334,
    day=totalday-304;
    Month=num2str(11);
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==334
    Month=num2str(11);
    Day=num2str(30);
end

if totalday>334 && totalday<365,
    day=totalday-334;
    Month=num2str(12);
    if day<10
        Day=strcat(num2str(0),num2str(day));
    else
        Day=num2str(day);
    end
end
if totalday==365
    Month=num2str(12);
    Day=num2str(31);
end

date=strcat(Month,'.',Day);


