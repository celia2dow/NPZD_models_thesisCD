% Plot mu_pT and beta_T against each other for a range of temperatures
close all
clear
temps = 10:0.5:30;
beta_20 = 0.66;
Q_10 = 1.88;
mu_20 = 0.03;
beta_T = beta_20 .* 1.066.^(temps-20);
mu_pT = mu_20 .*Q_10.^((temps-20)/10);

yyaxis left
plot(temps,beta_T)
ylabel("per day")
hold on
yyaxis right
plot(temps,mu_pT)
hold off
xlabel("Temperatures (^{\circ} C)")
ylabel("per day")
legend("\beta_T","\mu_{pT}")