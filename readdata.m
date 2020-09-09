load make.txt
load use.txt
disp (make);
disp (use);
coefmake = make./make(8,:);
coefuse = use./use(8,:);
disp (coefmake);
disp (coefuse);
coefmaketouse = make./use(8,:);
coefusetomake = use./make(8,:);