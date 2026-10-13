#include "poca.h"

#include "TPolyMarker3D.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TColor.h"
#include "TProfile3D.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TMarker.h"

POCACalculator::POCACalculator()
{
    ;
}

POCACalculator::POCACalculator(TString f_frecName, TString f_fpocaName)
{
    fRecName = f_frecName;
    fpocaName = f_fpocaName;
}

POCACalculator::POCACalculator(vector<TString> f_pocaslist, TString f_fpocaName)
{
    fpocalist = f_pocaslist;
    fpocaName = f_fpocaName;
}

POCACalculator::~POCACalculator()
{
    ;
}

bool POCACalculator::Init()
{
    for (int i = 0; i < 4; i++)
    {
        fRecFile[i] = new TFile(fpocalist[i], "read");
        if (!fRecFile[i]->IsOpen())
        {
            cout << " ERROR! " << fpocalist[i] << " cant open!" << endl;
            return false;
        }
        cout << fpocalist.size() << "," << fpocalist[3] << endl;
        fRecTree[i] = (TTree *)fRecFile[i]->Get("Track");
        fRecTree[i]->SetBranchAddress("event", &event);
        fRecTree[i]->SetBranchAddress("btofuseflag", &btofuseflag);
        fRecTree[i]->SetBranchAddress("X01fitPars", X01fitPars);
        fRecTree[i]->SetBranchAddress("Y01fitPars", Y01fitPars);
        fRecTree[i]->SetBranchAddress("X2345fitPars", X2345fitPars);
        fRecTree[i]->SetBranchAddress("Y2345fitPars", Y2345fitPars);
        fRecTree[i]->SetBranchAddress("useflag", &useflag);
        fRecTree[i]->SetBranchAddress("XfitPars", XfitPars);
        fRecTree[i]->SetBranchAddress("YfitPars", YfitPars);
    }

    fpocaFile = new TFile(fpocaName, "RECREATE");
    if (!fpocaFile->IsOpen())
    {
        cout << "ERROR! " << fpocaName << " cant open!" << endl;
        return false;
    }
    fpocaTree = new TTree("poca", "poca");
    fpocaTree->Branch("event", &event, "event/I");
    fpocaTree->Branch("btofuseflag", &btofuseflag, "btofuseflag/O");
    fpocaTree->Branch("X01fitPars", X01fitPars, "X01fitPars[3]/D");
    fpocaTree->Branch("Y01fitPars", Y01fitPars, "Y01fitPars[3]/D");
    fpocaTree->Branch("X2345fitPars", X2345fitPars, "X2345fitPars[3]/D");
    fpocaTree->Branch("Y2345fitPars", Y2345fitPars, "Y2345fitPars[3]/D");
    fpocaTree->Branch("pocax", &pocax, "pocax/D");
    fpocaTree->Branch("pocay", &pocay, "pocay/D");
    fpocaTree->Branch("pocaz", &pocaz, "pocaz/D");
    fpocaTree->Branch("scattering_angle", &scattering_angle, "scattering_angle/D");
    return true;
}

void POCACalculator::InitData()
{
    btofuseflag = 0;
    memset(X01fitPars, 0, sizeof(X01fitPars));
    memset(Y01fitPars, 0, sizeof(Y01fitPars));
    memset(X2345fitPars, 0, sizeof(X2345fitPars));
    memset(Y2345fitPars, 0, sizeof(Y2345fitPars));
    pocax = 0;
    pocay = 0;
    pocaz = 0;
    scattering_angle = 0;
}

void POCACalculator::Loop()
{

    TH2F *pocaxymapy;
    pocaxymapy = new TH2F(Form("pocazymapy"), Form("pocazymapy"), 1000, 3600, 4200, 100, -100, 100);
    pocaxymapy->GetXaxis()->SetTitle("Z [mm]");
    pocaxymapy->GetYaxis()->SetTitle("Y [mm]");

    TH2F *pocaxymapx;
    pocaxymapx = new TH2F(Form("pocazymapx"), Form("pocazymapx"), 1000, 3600, 4200, 100, -100, 100);
    pocaxymapx->GetXaxis()->SetTitle("Z [mm]");
    pocaxymapx->GetYaxis()->SetTitle("X [mm]");

    TH1F *angle;
    angle = new TH1F(Form("angle"), Form("angle"), 500, 0, 5);
    angle->GetXaxis()->SetTitle("angle");
    angle->GetYaxis()->SetTitle("counts");

    TH3F *voxel_hist = new TH3F("voxel_3d", "3D Voxel Density of POCA Points; X [mm]; Y [mm]; Z [mm]",
                                50, -100, 100,   // X轴: 50个bin, 3900-4100mm
                                50, -100, 100,   // Y轴: 50个bin, -100-100mm
                                50, 3900, 4100); // Z轴: 50个bin, 3900-4100mm

    // 创建散射角分布的二维数组
    std::vector<std::vector<TH1F *>> angle_distributions(25, std::vector<TH1F *>(25, nullptr));

    // 初始化所有分布直方图
    for (int i = 0; i < 25; i++)
    {
        for (int j = 0; j < 25; j++)
        {
            angle_distributions[i][j] = new TH1F(
                Form("angle_dist_z%d_x%d", i, j),
                Form("Scattering Angle Distribution (Z bin %d, X bin %d);Angle [deg];Counts", i, j),
                60, 0, 30);
        }
    }
    //////////////////////////
    const double z_min = 500.0, z_max = 9500.0; // Z范围 (mm)
    const double x_min = -100.0, x_max = 100.0; // X范围 (mm)
    const double y_min = -100.0, y_max = 100.0; // Y范围 (mm)
    const int nz_bins = 2000, nxy_bins = 200;   // bin数量，可调节分辨率

    TH2F *hXLines = new TH2F("hXLines", "Incident tracks in X-Z plane; Z (mm); X (mm)",
                             nz_bins, z_min, z_max, nxy_bins, x_min, x_max);
    TH2F *hYLines = new TH2F("hYLines", "Incident tracks in Y-Z plane; Z (mm); Y (mm)",
                             nz_bins, z_min, z_max, nxy_bins, y_min, y_max);

    int bin_num_xy = 40 / 5;
    int bin_num_z = 600 / 5;

    TH1F *sc_ang[bin_num_xy][bin_num_z];

    for (int i = 0; i < bin_num_xy; i++)
    {
        for (int j = 0; j < bin_num_z; j++)
        {
            sc_ang[i][j] = new TH1F(Form("angle%d_x%d", i, j), Form("angle%d_x%d", i, j), 50, -10, 10);
        }
    }

    TH2F *hSigma = new TH2F("hSigma",
                            "Gaussian #sigma of scattering angle vs (Z,Y); Z [mm]; Y [mm]",
                            bin_num_z, 3600, 4200, // Z bin 数、范围
                            bin_num_xy, -20, 20);  // Y bin 数、范围

    TH1F *y_noise = new TH1F("noise_y", "noise_y", 100, -50, 50);

    TH2F *h_xy = new TH2F("hxy",
                          "xy",
                          100, -50, 50,  // Z bin 数、范围
                          100, -50, 50); // Y bin 数、范围

    for (int kk = 0; kk < 4; kk++)
    {
        int nentries = fRecTree[kk]->GetEntriesFast();

        for (int i = 0; i < nentries; i++)
        {
            DrawProcessbar(i, nentries);
            InitData();
            fRecTree[kk]->GetEntry(i);
            if (btofuseflag)
            {
                TrackParams track_incoming(X01fitPars[0], X01fitPars[1], Y01fitPars[0], Y01fitPars[1]);
                TrackParams track_outgoing(X2345fitPars[0], X2345fitPars[1], Y2345fitPars[0], Y2345fitPars[1]);
                double angle_xz = TMath::ATan(track_outgoing.kx) - TMath::ATan(track_incoming.kx);
                double angle_yz = TMath::ATan(track_outgoing.ky) - TMath::ATan(track_incoming.ky);
                POCAResult result = CalculatePOCA(track_incoming, track_outgoing);
                pocax = result.poca_point.X();
                pocay = result.poca_point.Y();
                pocaz = result.poca_point.Z();
                scattering_angle = result.scattering_angle * 1000;
                // pocaxymapy->Fill(pocaz, pocay);
                // pocaxymapx->Fill(pocaz, pocax);
                if (TMath::Abs(angle_xz) * 1000 < 50 && TMath::Abs(angle_yz) * 1000 < 50)
                {
                    pocaxymapy->Fill(pocaz, pocay, TMath::Abs(angle_yz) * 100000);
                    pocaxymapx->Fill(pocaz, pocax, TMath::Abs(angle_xz) * 100000);
                    // pocaxymapy->Fill(pocaz, pocay);
                    // pocaxymapx->Fill(pocaz, pocax);
                }

                int binZ = pocaxymapy->GetXaxis()->FindBin(pocaz);
                int binX = pocaxymapy->GetYaxis()->FindBin(pocax);

                if (pocay < 20 && pocay > -20 && pocaz < 4200 && pocaz > 3600)
                {
                    int x_ind = pocay + 20;
                    int z_ind = pocaz - 3600;
                    x_ind = x_ind / 5;
                    z_ind = z_ind / 5;
                    sc_ang[x_ind][z_ind]->Fill(angle_yz * 1000);
                }

                if (pocaz > 3000 && pocaz < 3600)
                {
                    y_noise->Fill(pocay);
                }

                if (pocaz < 3972 && pocaz > 3964)
                {
                    h_xy->Fill(pocax, pocay);
                }

                // 填充对应的分布直方图
                if (binZ >= 1 && binZ <= 25 && binX >= 1 && binX <= 25)
                {
                    angle_distributions[binZ - 1][binX - 1]->Fill(scattering_angle);
                }
                voxel_hist->Fill(pocax, pocay, pocaz);
            }

            if (useflag)
            {
                double kx = XfitPars[0]; // dx/dz
                double bx = XfitPars[1]; // x at z=0
                double ky = YfitPars[0]; // dy/dz
                double by = YfitPars[1]; // y at z=0

                for (int iz = 1; iz <= nz_bins; iz++)
                {
                    double z = hXLines->GetXaxis()->GetBinCenter(iz);
                    double x = kx * z + bx;
                    double y = ky * z + by;

                    // 检查计算值是否在直方图范围内（防止填充溢出）
                    if (x >= x_min && x <= x_max)
                    {
                        hXLines->Fill(z, x);
                    }
                    if (y >= y_min && y <= y_max)
                    {
                        hYLines->Fill(z, y);
                    }
                }
            }

            fpocaTree->Fill();
        }
    }
    angle->Write();
    TCanvas *c3d = new TCanvas("c3d", "3D Voxel Density Plot", 1200, 800);
    voxel_hist->SetMinimum(10);
    voxel_hist->SetFillColorAlpha(kBlue, 0.1);
    voxel_hist->Draw("BOX2 Z");
    // voxel_hist->Draw("LEGO1");

    for (int i = 0; i < 25; i++)
    {
        for (int j = 0; j < 25; j++)
        {
            // if (angle_distributions[i][j]->GetEntries() > 100)
            // {
            //     // angle_distributions[i][j]->Fit("landau", "QR", "", 0, 15);
            //     angle_distributions[i][j]->Write();
            // }
            delete angle_distributions[i][j];
        }
    }

    for (int i = 0; i < bin_num_xy; i++)
    {
        for (int j = 0; j < bin_num_z; j++)
        {
            if (sc_ang[i][j]->GetEntries() > 200)
            {
                sc_ang[i][j]->Fit("gaus", "QR", "", -10, 10);
                double res = sc_ang[i][j]->GetFunction("gaus")->GetParameter(2);
                // sc_ang[i][j]->Write();

                double z_center = 3600 + j * 5 + 2.5; // bin 中心
                double y_center = -20 + i * 5 + 2.5;
                if (res * 1000 < 4000 && res * 1000 > 2000)
                    hSigma->Fill(z_center, y_center, res * 1000);
            }

            delete sc_ang[i][j];
        }
    }

    c3d->Update();
    c3d->Write();
    voxel_hist->Write();

    TCanvas *cSigma = new TCanvas("cSigma", "Sigma Map", 800, 600);
    cSigma->SetRightMargin(0.15);
    hSigma->SetStats(0);
    hSigma->SetZTitle("#sigma [deg]");
    hSigma->Draw("COLZ");
    cSigma->Write();
    hSigma->Write();

    TCanvas *cx = new TCanvas("cx", "x", 800, 600);
    pocaxymapx->Draw("COLZ");
    cx->Write();

    TCanvas *cy = new TCanvas("cy", "y", 800, 600);
    pocaxymapy->Draw("COLZ");
    cy->Write();

    TCanvas *ny = new TCanvas("ny", "ny", 800, 600);
    y_noise->Fit("gaus", "QR", "", -50, 50);
    y_noise->Draw();
    ny->Write();

    TCanvas *fhxy = new TCanvas("hxy", "hxf", 800, 600);
    h_xy->Draw();
    fhxy->Write();

    const int nbinsZ = pocaxymapy->GetNbinsX();
    const int nbinsY = pocaxymapy->GetNbinsY();

    // 存储每个 Z bin 的 Y 均值（索引 1..nbinsZ）
    vector<double> meanY(nbinsY + 1, 0.0);

    for (int iy = 1; iy <= nbinsY; ++iy)
    {
        double sumWY = 0.0, sumW = 0.0;
        for (int iz = 1; iz <= nbinsZ; ++iz)
        {
            double z_center = pocaxymapy->GetXaxis()->GetBinCenter(iz);
            if (z_center <= 3000.0 || z_center >= 3600.0)
                continue;
            double w = pocaxymapy->GetBinContent(iz, iy);
            if (w == 0.0)
                continue;
            sumWY += w;
            sumW++;
        }
        if (sumWY > 0.0)
            meanY[iy] = sumWY / sumW;
    }

    TH2F *pocaxymapy_sub = new TH2F(
        "pocaxymapy_sub",
        "POCA Y-Z (Y subtracted by Z-mean); Z [mm]; Y - <Y>(Z) [mm]",
        1000, 3600, 4200, // Z: 同原始
        100, -100, 100); // Y: bin 数增加/范围放宽

    // 第三步：将原始直方图内容移动到新直方图中（扣除 Y 均值）
    for (int iy = 1; iy <= nbinsY; ++iy)
    {
        double my = meanY[iy]; // 该 Z bin 的 Y 均值（3000-3600 外为 0）
        double z_center = pocaxymapy->GetXaxis()->GetBinCenter(iy);

        for (int iz = 1; iz <= nbinsZ; ++iz)
        {
            double content = pocaxymapy->GetBinContent(iz, iy);
            double error = pocaxymapy->GetBinError(iz, iy);
            double newcontent = content - my;
            if (newcontent <= 0)
                newcontent = 0;

            pocaxymapy_sub->SetBinContent(iz, iy, newcontent);

        }
    }

    // 第四步：保存或绘制新直方图
    TCanvas *cy_sub = new TCanvas("cy_sub", "Y subtracted", 800, 600);
    pocaxymapy_sub->Draw("COLZ");
    cy_sub->Write();
    pocaxymapy_sub->Write();

    delete pocaxymapy_sub;

    // --- 新增：处理并绘制X-Z和Y-Z直线密度图 ---
    TCanvas *cLines = new TCanvas("cLines", "Incident track projections", 1200, 500);
    cLines->Divide(2, 1);

    // X-Z平面
    cLines->cd(1);
    gPad->SetGrid();
    hXLines->SetStats(0);
    hXLines->Draw("COLZ");
    // 找出最大bin位置（即直线交点最密集的点）
    int maxBinX = hXLines->GetMaximumBin();
    int binXz, binXx, binXy;
    hXLines->GetBinXYZ(maxBinX, binXz, binXx, binXy);
    double z_peak_x = hXLines->GetXaxis()->GetBinCenter(binXz);
    double x_peak = hXLines->GetYaxis()->GetBinCenter(binXx);
    // 在图上标记该点
    TMarker *markX = new TMarker(z_peak_x, x_peak, 29); // 星形标记
    markX->SetMarkerColor(kRed);
    markX->SetMarkerSize(2);
    markX->Draw("SAME");
    // 添加文字
    TLatex *texX = new TLatex(z_peak_x, x_peak, Form("Peak: (%.1f, %.1f)", z_peak_x, x_peak));
    texX->SetTextColor(kRed);
    texX->SetTextSize(0.03);
    texX->Draw("SAME");

    // Y-Z平面
    cLines->cd(2);
    gPad->SetGrid();
    hYLines->SetStats(0);
    hYLines->Draw("COLZ");
    int maxBinY = hYLines->GetMaximumBin();
    int binYz, binYy, binYtmp;
    hYLines->GetBinXYZ(maxBinY, binYz, binYy, binYtmp);
    double z_peak_y = hYLines->GetXaxis()->GetBinCenter(binYz);
    double y_peak = hYLines->GetYaxis()->GetBinCenter(binYy);
    TMarker *markY = new TMarker(z_peak_y, y_peak, 29);
    markY->SetMarkerColor(kRed);
    markY->SetMarkerSize(2);
    markY->Draw("SAME");
    TLatex *texY = new TLatex(z_peak_y, y_peak, Form("Peak: (%.1f, %.1f)", z_peak_y, y_peak));
    texY->SetTextColor(kRed);
    texY->SetTextSize(0.03);
    texY->Draw("SAME");

    cLines->Write();  // 保存画布
    hXLines->Write(); // 保存直方图
    hYLines->Write();

    delete voxel_hist;
    delete c3d;
    delete pocaxymapx;
    delete pocaxymapy;
    delete angle;
    delete cLines;
    delete hXLines;
    delete hYLines;
    fpocaTree->Write();
    fpocaFile->Flush();
    fpocaFile->Close();
}

POCACalculator::POCAResult POCACalculator::CalculatePOCA(const TrackParams &track1, const TrackParams &track2)
{
    POCAResult result;
    double kx1 = track1.kx, bx1 = track1.bx, ky1 = track1.ky, by1 = track1.by;
    double kx2 = track2.kx, bx2 = track2.bx, ky2 = track2.ky, by2 = track2.by;
    TMatrixD A(2, 2);
    TVectorD B(2);
    A(0, 0) = kx1 * kx1 + ky1 * ky1 + 1;
    A(0, 1) = -(kx1 * kx2 + ky1 * ky2 + 1);
    A(1, 0) = -(kx1 * kx2 + ky1 * ky2 + 1);
    A(1, 1) = kx2 * kx2 + ky2 * ky2 + 1;
    B(0) = kx1 * (bx2 - bx1) + ky1 * (by2 - by1);
    B(1) = kx2 * (bx1 - bx2) + ky2 * (by1 - by2);
    try
    {
        TMatrixD A_inv = A;
        A_inv.Invert();
        TVectorD solution = A_inv * B;
        double z1 = solution(0);
        double z2 = solution(1);
        result.point_on_track1.SetXYZ(kx1 * z1 + bx1, ky1 * z1 + by1, z1);
        result.point_on_track2.SetXYZ(kx2 * z2 + bx2, ky2 * z2 + by2, z2);
        result.poca_point = (result.point_on_track1 + result.point_on_track2) * 0.5;
        result.distance = (result.point_on_track1 - result.point_on_track2).Mag();
        result.scattering_angle = CalculateScatteringAngle(track1, track2);
    }
    catch (...)
    {
        std::cout << "Warning: Matrix inversion failed, using fallback method for parallel tracks." << std::endl;
        return HandleParallelTracks(track1, track2);
    }
    return result;
}

double POCACalculator::CalculateScatteringAngle(const TrackParams &track1, const TrackParams &track2)
{
    TVector3 dir1(track1.kx, track1.ky, 1);
    TVector3 dir2(track2.kx, track2.ky, 1);

    dir1 = dir1.Unit();
    dir2 = dir2.Unit();

    double dot_product = dir1.Dot(dir2);
    // 防止浮点误差导致数值问题
    dot_product = TMath::Max(-1.0, TMath::Min(1.0, dot_product));

    return TMath::ACos(dot_product);
}

void POCACalculator::CalculateScatteringComponents(const TrackParams &track1, const TrackParams &track2, double &angle_x, double &angle_y)
{
    // X-Z 平面上的角度
    double theta1_x = TMath::ATan(track1.kx);
    double theta2_x = TMath::ATan(track2.kx);
    angle_x = theta2_x - theta1_x;

    // Y-Z 平面上的角度
    double theta1_y = TMath::ATan(track1.ky);
    double theta2_y = TMath::ATan(track2.ky);
    angle_y = theta2_y - theta1_y;
}

POCACalculator::POCAResult POCACalculator::HandleParallelTracks(const TrackParams &track1, const TrackParams &track2)
{
    POCAResult result;

    // 选择一个固定的 z 值（例如 z=0）
    double z_fixed = 0;

    result.point_on_track1.SetXYZ(
        track1.kx * z_fixed + track1.bx,
        track1.ky * z_fixed + track1.by,
        z_fixed);

    result.point_on_track2.SetXYZ(
        track2.kx * z_fixed + track2.bx,
        track2.ky * z_fixed + track2.by,
        z_fixed);

    result.poca_point = (result.point_on_track1 + result.point_on_track2) * 0.5;
    result.distance = (result.point_on_track1 - result.point_on_track2).Mag();
    result.scattering_angle = CalculateScatteringAngle(track1, track2);

    return result;
}