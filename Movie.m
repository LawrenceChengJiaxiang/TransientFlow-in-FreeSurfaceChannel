% 初始化一个电影矩阵
M = moviein(432);
Qin = [];
% 创建电影
for k = 1:432
   Qin = [Qin, Qinput(k*12/432)];
   subplot(2,1,1);plot(0:12/432:12*(k-1)/432,Qin(1:k));
   xlim([0 12]);ylim([0 1])
   subplot(2,1,2);plot(0:100:7000,Hmovie(k*100+1,1:71));
   xlim([0 7000]);ylim([0.25 1])
    % 调用getframe函数生成每个帧
   M(k) = getframe;
end
% 调用movie函数将电影动画矩阵M(k)播放5次
movie(M,1);
