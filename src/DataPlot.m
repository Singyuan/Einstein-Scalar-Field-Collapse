function data = DataPlot(data, recorder, totaltime, alphaAdap)
data.centralalpha(recorder) = alphaAdap(1);
data.timecoord(recorder) = totaltime;
end