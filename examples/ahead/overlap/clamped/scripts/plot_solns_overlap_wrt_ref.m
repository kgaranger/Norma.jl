close all; clear all;

coords1 = csvread('clamped-1-refe.csv');
N1 = length(coords1);
n1 = N1/3;
x1 = coords1(:,1);
y1 = coords1(:,2);
z1 = coords1(:,3);
ind1 = find(x1 == min(x1) & y1 == min(y1));
coords2 = csvread('clamped-2-refe.csv');
N2 = length(coords2);
n2 = N2/3;
x2 = coords2(:,1);
y2 = coords2(:,2);
z2 = coords2(:,3);
ind2 = find(x2 == min(x2) & y2 == min(y2));
z = sort(unique([z1;z2]));
dispz1 = []; veloz1 = []; accez1 = [];
dispz2 = []; veloz2 = []; accez2 = [];
dispzRef1 = []; velozRef1 = []; accezRef1 = [];
dispzRef2 = []; velozRef2 = []; accezRef2 = [];
disp_computed = []; velo_computed = []; acce_computed = [];
disp_exact = []; velo_exact = []; acce_exact = [];
fig = figure();
save_figs = 1;
ctr = 1;
for i=0:100:10000
  if (i < 10)
    disp1_file_name = strcat('clamped-1-disp-000', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-000', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-000', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-000', num2str(i), '.csv');
    disp2_file_name = strcat('clamped-2-disp-000', num2str(i), '.csv');
    velo2_file_name = strcat('clamped-2-velo-000', num2str(i), '.csv');
    acce2_file_name = strcat('clamped-2-acce-000', num2str(i), '.csv');
    time2_file_name = strcat('clamped-2-time-000', num2str(i), '.csv');
    dispRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-000', num2str(i), '.csv');
    veloRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-000', num2str(i), '.csv');
    acceRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-000', num2str(i), '.csv');
    dispRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-disp-000', num2str(i), '.csv');
    veloRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-velo-000', num2str(i), '.csv');
    acceRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-acce-000', num2str(i), '.csv');
  elseif (i < 100)
    disp1_file_name = strcat('clamped-1-disp-00', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-00', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-00', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-00', num2str(i), '.csv');
    disp2_file_name = strcat('clamped-2-disp-00', num2str(i), '.csv');
    velo2_file_name = strcat('clamped-2-velo-00', num2str(i), '.csv');
    acce2_file_name = strcat('clamped-2-acce-00', num2str(i), '.csv');
    time2_file_name = strcat('clamped-2-time-00', num2str(i), '.csv');
    dispRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-00', num2str(i), '.csv');
    veloRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-00', num2str(i), '.csv');
    acceRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-00', num2str(i), '.csv');
    dispRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-disp-00', num2str(i), '.csv');
    veloRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-velo-00', num2str(i), '.csv');
    acceRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-acce-00', num2str(i), '.csv');
  elseif (i < 1000)
    disp1_file_name = strcat('clamped-1-disp-0', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-0', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-0', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-0', num2str(i), '.csv');
    disp2_file_name = strcat('clamped-2-disp-0', num2str(i), '.csv');
    velo2_file_name = strcat('clamped-2-velo-0', num2str(i), '.csv');
    acce2_file_name = strcat('clamped-2-acce-0', num2str(i), '.csv');
    time2_file_name = strcat('clamped-2-time-0', num2str(i), '.csv');
    dispRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-0', num2str(i), '.csv');
    veloRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-0', num2str(i), '.csv');
    acceRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-0', num2str(i), '.csv');
    dispRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-disp-0', num2str(i), '.csv');
    veloRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-velo-0', num2str(i), '.csv');
    acceRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-acce-0', num2str(i), '.csv');
  else
    disp1_file_name = strcat('clamped-1-disp-', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-', num2str(i), '.csv');
    disp2_file_name = strcat('clamped-2-disp-', num2str(i), '.csv');
    velo2_file_name = strcat('clamped-2-velo-', num2str(i), '.csv');
    acce2_file_name = strcat('clamped-2-acce-', num2str(i), '.csv');
    time2_file_name = strcat('clamped-2-time-', num2str(i), '.csv');
    dispRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-', num2str(i), '.csv');
    veloRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-', num2str(i), '.csv');
    acceRef1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-', num2str(i), '.csv');
    dispRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-disp-', num2str(i), '.csv');
    veloRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-velo-', num2str(i), '.csv');
    acceRef2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-2-acce-', num2str(i), '.csv');
  end
  d1 = dlmread(disp1_file_name);
  v1 = dlmread(velo1_file_name);
  a1 = dlmread(acce1_file_name);
  t1 = dlmread(time1_file_name);
  dispz1 = [dispz1, d1(:,3)];
  veloz1 = [veloz1, v1(:,3)];
  accez1 = [accez1, a1(:,3)];
  d2 = dlmread(disp2_file_name);
  v2 = dlmread(velo2_file_name);
  a2 = dlmread(acce2_file_name);
  t2 = dlmread(time2_file_name);
  dispz2 = [dispz2, d2(:,3)];
  veloz2 = [veloz2, v2(:,3)];
  accez2 = [accez2, a2(:,3)];
  dRef1 = dlmread(dispRef1_file_name);
  vRef1 = dlmread(veloRef1_file_name);
  aRef1 = dlmread(acceRef1_file_name);
  dispzRef1 = [dispzRef1, dRef1(:,3)];
  velozRef1 = [velozRef1, vRef1(:,3)];
  accezRef1 = [accezRef1, aRef1(:,3)];
  dRef2 = dlmread(dispRef2_file_name);
  vRef2 = dlmread(veloRef2_file_name);
  aRef2 = dlmread(acceRef2_file_name);
  dispzRef2 = [dispzRef2, dRef2(:,3)];
  velozRef2 = [velozRef2, vRef2(:,3)];
  accezRef2 = [accezRef2, aRef2(:,3)];
  z1ind1 = z1(ind1);
  z1ind1 = z1(ind1);
  z2ind2 = z2(ind2);
  zz = uniquetol(sort([z1ind1;z2ind2]), 1e-5);
  dz1ind1 = dispz1(ind1,ctr);
  dz2ind2 = dispz2(ind2,ctr);
  vz1ind1 = veloz1(ind1,ctr);
  vz2ind2 = veloz2(ind2,ctr);
  az1ind1 = accez1(ind1,ctr);
  az2ind2 = accez2(ind2,ctr);
  dzRef1ind1 = dispzRef1(ind1,ctr);
  dzRef2ind2 = dispzRef2(ind2,ctr);
  vzRef1ind1 = velozRef1(ind1,ctr);
  vzRef2ind2 = velozRef2(ind2,ctr);
  azRef1ind1 = accezRef1(ind1,ctr);
  azRef2ind2 = accezRef2(ind2,ctr);
  for k = 1:length(zz)
    ii1 = find(z1ind1 == zz(k));
    ii2 = find(z2ind2 == zz(k));
    dispz_merged(k,1) = 0;
    veloz_merged(k,1) = 0;
    accez_merged(k,1) = 0;
    dispzRef_merged(k,1) = 0;
    velozRef_merged(k,1) = 0;
    accezRef_merged(k,1) = 0;
    if length(ii1) > 0
      dispz_merged(k,1) = dispz_merged(k,1) + dz1ind1(ii1);
      veloz_merged(k,1) = veloz_merged(k,1) + vz1ind1(ii1);
      accez_merged(k,1) = accez_merged(k,1) + az1ind1(ii1);
      dispzRef_merged(k,1) = dispzRef_merged(k,1) + dzRef1ind1(ii1);
      velozRef_merged(k,1) = velozRef_merged(k,1) + vzRef1ind1(ii1);
      accezRef_merged(k,1) = accezRef_merged(k,1) + azRef1ind1(ii1);
    end
    if length(ii2) > 0
      dispz_merged(k,1) = dispz_merged(k,1) + dz2ind2(ii2);
      veloz_merged(k,1) = veloz_merged(k,1) + vz2ind2(ii2);
      accez_merged(k,1) = accez_merged(k,1) + az2ind2(ii2);
      dispzRef_merged(k,1) = dispzRef_merged(k,1) + dzRef2ind2(ii2);
      velozRef_merged(k,1) = velozRef_merged(k,1) + vzRef2ind2(ii2);
      accezRef_merged(k,1) = accezRef_merged(k,1) + azRef2ind2(ii2);
    end
    if length(ii1) + length(ii2) > 1
      dispz_merged(k,1) = dispz_merged(k,1)/2;
      veloz_merged(k,1) = veloz_merged(k,1)/2;
      accez_merged(k,1) = accez_merged(k,1)/2;
      dispzRef_merged(k,1) = dispzRef_merged(k,1)/2;
      velozRef_merged(k,1) = velozRef_merged(k,1)/2;
      accezRef_merged(k,1) = accezRef_merged(k,1)/2;
    end
  end
  disp_computed = [disp_computed, dispz_merged];
  velo_computed = [velo_computed, veloz_merged];
  acce_computed = [acce_computed, accez_merged];
  disp_exact = [disp_exact, dispzRef_merged];
  velo_exact = [velo_exact, velozRef_merged];
  acce_exact = [acce_exact, accezRef_merged];

  subplot(3,1,1);
  ax = gca;
  plot(z1(ind1), dispz1(ind1,ctr), '-b');
  hold on;
  plot(z2(ind2), dispz2(ind2,ctr), '-r');
  hold on;
  plot(z1(ind1), dispzRef1(ind1,ctr), '--c');
  hold on;
  plot(z2(ind2), dispzRef2(ind2,ctr), '--g');
  xlabel('z');
  ylabel('z-disp');
  hold on;
  title(['displacement snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -0.001 0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,2);
  ax = gca;
  plot(z1(ind1), veloz1(ind1,ctr), '-b');
  hold on;
  plot(z2(ind2), veloz2(ind2,ctr), '-r');
  hold on;
  plot(z1(ind1), velozRef1(ind1,ctr), '--c');
  hold on;
  plot(z2(ind2), velozRef2(ind2,ctr), '--g');
  xlabel('z');
  ylabel('z-velo');
  hold on;
  title(['velocity snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -3e4*0.001 3e4*0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,3);
  ax = gca;
  plot(z1(ind1), accez1(ind1,ctr), '-b');
  hold on;
  plot(z2(ind2), accez2(ind2,ctr), '-r');
  hold on;
  plot(z1(ind1), accezRef1(ind1,ctr), '--c');
  hold on;
  plot(z2(ind2), accezRef2(ind2,ctr), '--g');
  xlabel('z');
  ylabel('z-acce');
  title(['acceleration snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -2.5e9*0.001 2.5e9*0.001]);
  hold on;
  ax.NextPlot = 'replaceChildren';
  pause(0.5)
  %pause()
  if (save_figs == 1)
    if (ctr < 10)
      filename = strcat('soln_000', num2str(ctr), '.png');
      filename2 = strcat('soln_000', num2str(ctr), '.fig');
    elseif (ctr < 100)
      filename = strcat('soln_00', num2str(ctr), '.png');
      filename2 = strcat('soln_00', num2str(ctr), '.fig');
    elseif (ctr < 1000)
      filename = strcat('soln_0', num2str(ctr), '.png');
      filename2 = strcat('soln_0', num2str(ctr), '.fig');
    else
      filename = strcat('soln_', num2str(ctr), '.png');
      filename2 = strcat('soln_', num2str(ctr), '.fig');
    end
    saveas(fig,filename)
    saveas(fig,filename2)
  end
  ctr = ctr + 1;
end

sz = size(disp_exact);
numerator_disp = 0;
denomenator_disp = 0;
numerator_velo = 0;
denomenator_velo = 0;
numerator_acce = 0;
denomenator_acce = 0;
for i=1:sz(2)
  numerator_disp = numerator_disp + norm(disp_computed(:,i) - disp_exact(:,i))^2;
  denomenator_disp = denomenator_disp + norm(disp_exact(:,i))^2;
  numerator_velo = numerator_velo + norm(velo_computed(:,i) - velo_exact(:,i))^2;
  denomenator_velo = denomenator_velo + norm(velo_exact(:,i))^2;
  numerator_acce = numerator_acce + norm(acce_computed(:,i) - acce_exact(:,i))^2;
  denomenator_acce = denomenator_acce + norm(acce_exact(:,i))^2;
end
dispz_relerr = sqrt(numerator_disp/denomenator_disp);
veloz_relerr = sqrt(numerator_velo/denomenator_velo);
accez_relerr = sqrt(numerator_acce/denomenator_acce);

fprintf('z-disp avg rel error = %e\n', dispz_relerr);
fprintf('z-velo avg rel error = %e\n', veloz_relerr);
fprintf('z-acce avg rel error = %e\n', accez_relerr);
