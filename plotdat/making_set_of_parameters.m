load ../result/parameter5_step0.dat;
load ../result/parameter5_step1.dat;
load ../result/parameter5_step2.dat;
load ../result/parameter5_step3.dat;
load ../result/parameter5_step4.dat;
load ../result/parameter5_step5.dat;

load ../result/parameter6_step0.dat;
load ../result/parameter6_step1.dat;
load ../result/parameter6_step2.dat;
load ../result/parameter6_step3.dat;
load ../result/parameter6_step4.dat;
load ../result/parameter6_step5.dat;

load ../result/parameter7_step0.dat;
load ../result/parameter7_step1.dat;
load ../result/parameter7_step2.dat;
load ../result/parameter7_step3.dat;
load ../result/parameter7_step4.dat;
load ../result/parameter7_step5.dat;

x = parameter5_step0(:,1);
xs =[0:1:350];

yp5s0 = parameter5_step0(:,2);
yp5s1 = parameter5_step1(:,2);
yp5s2 = parameter5_step2(:,2);
yp5s3 = parameter5_step3(:,2);
yp5s4 = parameter5_step4(:,2);
yp5s5 = parameter5_step5(:,2);

yp6s0 = parameter6_step0(:,2);
yp6s1 = parameter6_step1(:,2);
yp6s2 = parameter6_step2(:,2);
yp6s3 = parameter6_step3(:,2);
yp6s4 = parameter6_step4(:,2);
yp6s5 = parameter6_step5(:,2);

yp7s0 = parameter7_step0(:,2);
yp7s1 = parameter7_step1(:,2);
yp7s2 = parameter7_step2(:,2);
yp7s3 = parameter7_step3(:,2);
yp7s4 = parameter7_step4(:,2);
yp7s5 = parameter7_step5(:,2);

yintp5s0 = interp1(x, yp5s0, xs, "pchip");
yintp5s1 = interp1(x, yp5s1, xs, "pchip");
yintp5s2 = interp1(x, yp5s2, xs, "pchip");
yintp5s3 = interp1(x, yp5s3, xs, "pchip");
yintp5s4 = interp1(x, yp5s4, xs, "pchip");
yintp5s5 = interp1(x, yp5s5, xs, "pchip");

yintp6s0 = interp1(x, yp6s0, xs, "pchip");
yintp6s1 = interp1(x, yp6s1, xs, "pchip");
yintp6s2 = interp1(x, yp6s2, xs, "pchip");
yintp6s3 = interp1(x, yp6s3, xs, "pchip");
yintp6s4 = interp1(x, yp6s4, xs, "pchip");
yintp6s5 = interp1(x, yp6s5, xs, "pchip");

yintp7s0 = interp1(x, yp7s0, xs, "pchip");
yintp7s1 = interp1(x, yp7s1, xs, "pchip");
yintp7s2 = interp1(x, yp7s2, xs, "pchip");
yintp7s3 = interp1(x, yp7s3, xs, "pchip");
yintp7s4 = interp1(x, yp7s4, xs, "pchip");
yintp7s5 = interp1(x, yp7s5, xs, "pchip");

p5s0 = [xs; yintp5s0];
p5s1 = [xs; yintp5s1];
p5s2 = [xs; yintp5s2];
p5s3 = [xs; yintp5s3];
p5s4 = [xs; yintp5s4];
p5s5 = [xs; yintp5s5];

p6s0 = [xs; yintp6s0];
p6s1 = [xs; yintp6s1];
p6s2 = [xs; yintp6s2];
p6s3 = [xs; yintp6s3];
p6s4 = [xs; yintp6s4];
p6s5 = [xs; yintp6s5];

p7s0 = [xs; yintp7s0];
p7s1 = [xs; yintp7s1];
p7s2 = [xs; yintp7s2];
p7s3 = [xs; yintp7s3];
p7s4 = [xs; yintp7s4];
p7s5 = [xs; yintp7s5];

p5s0 = p5s0';
p5s1 = p5s1';
p5s2 = p5s2';
p5s3 = p5s3';
p5s4 = p5s4';
p5s5 = p5s5';

p6s0 = p6s0';
p6s1 = p6s1';
p6s2 = p6s2';
p6s3 = p6s3';
p6s4 = p6s4';
p6s5 = p6s5';

p7s0 = p7s0';
p7s1 = p7s1';
p7s2 = p7s2';
p7s3 = p7s3';
p7s4 = p7s4';
p7s5 = p7s5';

save p5s0.dat p5s0;
save p5s1.dat p5s1;
save p5s2.dat p5s2;
save p5s3.dat p5s3;
save p5s4.dat p5s4;
save p5s5.dat p5s5;

save p6s0.dat p6s0;
save p6s1.dat p6s1;
save p6s2.dat p6s2;
save p6s3.dat p6s3;
save p6s4.dat p6s4;
save p6s5.dat p6s5;

save p7s0.dat p7s0;
save p7s1.dat p7s1;
save p7s2.dat p7s2;
save p7s3.dat p7s3;
save p7s4.dat p7s4;
save p7s5.dat p7s5;

cp5s0 = [255 200 200] ./ 255;
cp5s1 = [255 160 160] ./ 255;
cp5s2 = [255 120 120] ./ 255;
cp5s3 = [255 80 80] ./ 255;
cp5s4 = [255 40 40] ./ 255;
cp5s5 = [255 0 0] ./ 255;

cp6s0 = [200 255 200] ./ 255;
cp6s1 = [160 255 160] ./ 255;
cp6s2 = [120 255 120] ./ 255;
cp6s3 = [80 255 80] ./ 255;
cp6s4 = [40 255 40] ./ 255;
cp6s5 = [0 255 0] ./ 255;

cp7s0 = [200 200 255] ./ 255;
cp7s1 = [160 160 255] ./ 255;
cp7s2 = [120 120 255] ./ 255;
cp7s3 = [80 80 255] ./ 255;
cp7s4 = [40 40 255] ./ 255;
cp7s5 = [0 0 255] ./ 255;

set (gcf, "papersize", [3.5, 3.5]) 
saturation = figure('Position',[0,0,1800,1800]);

subplot(3, 1, 1);
plot(xs, yintp5s0, "color", cp5s0, "linewidth", 17,
  xs, yintp5s1, "color", cp5s1, "linewidth", 14,
  xs, yintp5s2, "color", cp5s2, "linewidth", 11,
  xs, yintp5s3, "color", cp5s3, "linewidth", 8,
  xs, yintp5s4, "color", cp5s4, "linewidth", 5,
  xs, yintp5s5, "color", cp5s5, "linewidth", 2);
axes=get(gcf, "currentaxes");
set(axes, "fontsize", 20, "linewidth", 2);
axis([0, 350, 0, 1]);
legend_parameters = legend ("0 день", "73 день", "146 день", "219 день", "365 день");
legend (legend_parameters, "location", "northeastoutside");
set (legend_parameters, "fontsize", 30);

subplot(3, 1, 2);
plot(xs, yintp6s0, "color", cp6s0, "linewidth", 17,
  xs, yintp6s1, "color", cp6s1, "linewidth", 14,
  xs, yintp6s2, "color", cp6s2, "linewidth", 11,
  xs, yintp6s3, "color", cp6s3, "linewidth", 8,
  xs, yintp6s4, "color", cp6s4, "linewidth", 5,
  xs, yintp6s5, "color", cp6s5, "linewidth", 2);
axes=get(gcf, "currentaxes");
set(axes, "fontsize", 20, "linewidth", 2);
legend_parameters = legend ("0 день", "73 день", "146 день", "219 день", "365 день");
legend (legend_parameters, "location", "northeastoutside");
set (legend_parameters, "fontsize", 30);

subplot(3, 1, 3);
plot(xs, yintp7s0, "color", cp7s0, "linewidth", 17,
  xs, yintp7s1, "color", cp7s1, "linewidth", 14,
  xs, yintp7s2, "color", cp7s2, "linewidth", 11,
  xs, yintp7s3, "color", cp7s3, "linewidth", 8,
  xs, yintp7s4, "color", cp7s4, "linewidth", 5,
  xs, yintp7s5, "color", cp7s5, "linewidth", 2);
axes=get(gcf, "currentaxes");
set(axes, "fontsize", 20, "linewidth", 2);
legend_parameters = legend ("0 день", "73 день", "146 день", "219 день", "365 день");
legend (legend_parameters, "location", "northeastoutside");
set (legend_parameters, "fontsize", 30);

saveas(1, "saturation.png");