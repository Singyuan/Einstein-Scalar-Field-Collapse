function VideoMaker(rDense, sfcap, totaltime)
plot([-fliplr(rDense(2:end)) rDense], [fliplr(sfcap(1, 2:end)) sfcap(1, :)],'LineWidth',1.5)
xlabel('r')
axis([-inf, inf, -0.15, 0.15])
tit = ['time = ', num2str(totaltime)];
title(tit)
drawnow
end