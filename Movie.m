% ��ʼ��һ����Ӱ����
M = moviein(432);
Qin = [];
% ������Ӱ
for k = 1:432
   Qin = [Qin, Qinput(k*12/432)];
   subplot(2,1,1);plot(0:12/432:12*(k-1)/432,Qin(1:k));
   xlim([0 12]);ylim([0 1])
   subplot(2,1,2);plot(0:100:7000,Hmovie(k*100+1,1:71));
   xlim([0 7000]);ylim([0.25 1])
    % ����getframe��������ÿ��֡
   M(k) = getframe;
end
% ����movie��������Ӱ��������M(k)����5��
movie(M,1);
