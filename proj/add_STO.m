function y_STO=add_STO(y,iSTO)
    if iSTO>=0
        y_STO = [y(iSTO+1:end) zeros(1,iSTO)]; % advance|��ǰ
    else
        y_STO = [zeros(1,-iSTO) y(1:end+iSTO)];  % delay|�ͺ�
    end