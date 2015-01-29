%%Bu fonksiyon, GROUND TRUTH ile elde edilen sonuc görüntüsü için Precision
%%ve Recall degerlerini hesaplar. Ayrica GOSTER girisiyle de sonucu renkli
%%bir sekilde kullaniciya gösterir.
%v2'den farki pre ve rec'den farkli diger performans kriterlerinin de
%hesaba katilmasidir.

function [Precision Recall QP Specificity Accuracy Fmeasure] = fPreRecallHesapla_v3(GT,SONUC,GOSTER,FIG_NO)


%GT = logical(GT);
%SONUC = logical(SONUC);


% TPM = GT & SONUC;           %True Positive Matrix
% FPM = not(GT) & SONUC;      %False Positive Matrix
% FNM = GT & not(SONUC);      %False Negative Matrix
% TNM = not(GT) & not(SONUC);

TPM = GT .* SONUC;%True Positive Matrix
FPM = (1 - GT) .* SONUC;%False Positive Matrix
FNM = GT .* (1 - SONUC);%False Negative Matrix
TNM = (1 - GT) .* (1 - SONUC);

% figure(191),imshow(TPM,[])
%  figure(192),imshow(FPM,[])
%  figure(193),imshow(FNM,[])

TruePositive = sum(single(TPM(:)));
FalsePositive = sum(single(FPM(:)));
FalseNegative = sum(single(FNM(:)));
TrueNegative = sum(single(TNM(:)));
TP = TruePositive;
TN = TrueNegative;
FP = FalsePositive;
FN = FalseNegative;

Precision = TP /(TP + FP);%Correctness
Recall = TP / (TP + FN);  %Completeness
QP = TP /(TP + FN + FP);  %Quality Present?
Specificity = TN/(TN + FP);   %True Negative Rate
Accuracy = (TP + TN)/(TP + TN + FN + FP);%bence en iyi olcum bu degil...
Fmeasure = 2*Precision*Recall/(Precision+Recall);%Harmonic Mean of Precision and Recall 

if isnan(Precision),   Precision = 0;   end
if isnan(Recall),      Recall = 0;      end
if isnan(QP),          QP = 0;          end
if isnan(Accuracy),    Accuracy = 0;    end
if isnan(Specificity), Specificity = 0; end
if isnan(Fmeasure),    Fmeasure = 0;    end

if strcmpi(GOSTER,'goster')
  
    GORUNTU = single(zeros(size(GT)));
    GORUNTU(:,:,1) = FNM;
    GORUNTU(:,:,2) = TPM;
    GORUNTU(:,:,3) = FPM;
    figure(FIG_NO),imshow(GORUNTU,[]);
%     title(sprintf(['Pre = ' num2str(Precision) ' - Rec = ' num2str(Recall) ...
%                    'Q.F.= ' num2str(QP) ' - Acc = ' num2str(Accuracy) ...
%                    'Spe = ' num2str(Specificity) ' - Fmeas = ' num2str(Fmeasure) ]));
title(sprintf('Pre = %2.3g, Rec = %2.3g, Q_f = %2.3g, Fmeas = %2.3g, Acc = %2.3g, Spe = %2.3g',Precision,Recall,QP,Fmeasure,Accuracy,Specificity));
end
