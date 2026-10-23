# ============================================================
# Script PowerShell - Lancement des expériences PageRank
# ============================================================
# Usage : .\run_experiments.ps1
# Adapter les chemins $ExePath et $MatricesDir selon votre config
# ============================================================

# --- CONFIGURATION ---
$ExePath    = ".\main.exe"          # Chemin vers l'exécutable
$MatricesDir = "."                       # Dossier contenant les .mtx
$Matrix1    = "Harvard500.mtx"
$Matrix2    = "web-Google.mtx"
$Matrix3    = "uk-2002.mtx"
$LogFile    = "run_log.txt"

# --- OPTIMISATIONS SYSTÈME ---
Write-Host "=== Configuration du système ===" -ForegroundColor Cyan

# Mode haute performance
powercfg /setactive SCHEME_MIN
Write-Host "[OK] Mode haute performance activé"

# Désactiver la mise en veille
powercfg /change standby-timeout-ac 0
powercfg /change monitor-timeout-ac 0
Write-Host "[OK] Mise en veille désactivée"

# Priorité haute pour le processus
$CurrentProcess = Get-Process -Id $PID
$CurrentProcess.PriorityClass = "High"
Write-Host "[OK] Priorité haute définie"

# Afficher les infos système dans le log
"=== Expériences PageRank ===" | Out-File $LogFile
"Date de lancement : $(Get-Date)" | Out-File $LogFile -Append
"Machine : $env:COMPUTERNAME" | Out-File $LogFile -Append
"OS : $((Get-WmiObject Win32_OperatingSystem).Caption)" | Out-File $LogFile -Append
"CPU : $((Get-WmiObject Win32_Processor).Name)" | Out-File $LogFile -Append
"RAM : $([math]::Round((Get-WmiObject Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 2)) GB" | Out-File $LogFile -Append
"" | Out-File $LogFile -Append

Write-Host ""
Write-Host "=== Infos machine ===" -ForegroundColor Cyan
Write-Host "CPU : $((Get-WmiObject Win32_Processor).Name)"
Write-Host "RAM : $([math]::Round((Get-WmiObject Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 2)) GB"
Write-Host "OS  : $((Get-WmiObject Win32_OperatingSystem).Caption)"
Write-Host ""

# --- LANCEMENT HARVARD500 ---
Write-Host "=== Lancement Harvard500 ===" -ForegroundColor Green
"=== Harvard500 ===" | Out-File $LogFile -Append
$Start1 = Get-Date
"Début : $Start1" | Out-File $LogFile -Append

& $ExePath "$MatricesDir\$Matrix1"

$End1 = Get-Date
$Duration1 = $End1 - $Start1
"Fin : $End1" | Out-File $LogFile -Append
"Durée : $($Duration1.ToString())" | Out-File $LogFile -Append
"" | Out-File $LogFile -Append

Write-Host "Harvard500 terminé en $($Duration1.TotalMinutes.ToString('F1')) minutes" -ForegroundColor Green
Write-Host ""

# --- LANCEMENT WEB-GOOGLE ---
Write-Host "=== Lancement web-Google ===" -ForegroundColor Green
"=== web-Google ===" | Out-File $LogFile -Append
$StartG = Get-Date
"Début : $StartG" | Out-File $LogFile -Append

& $ExePath "$MatricesDir\$Matrix2"

$EndG = Get-Date
$DurationG = $EndG - $StartG
"Fin : $EndG" | Out-File $LogFile -Append
"Durée : $($DurationG.ToString())" | Out-File $LogFile -Append
"" | Out-File $LogFile -Append

Write-Host "web-Google terminé en $($DurationG.TotalMinutes.ToString('F1')) minutes" -ForegroundColor Green
Write-Host ""

# --- LANCEMENT UK-2002 ---
Write-Host "=== Lancement uk-2002 ===" -ForegroundColor Yellow
"=== uk-2002 ===" | Out-File $LogFile -Append
$Start2 = Get-Date
"Début : $Start2" | Out-File $LogFile -Append

& $ExePath "$MatricesDir\$Matrix3"

$End2 = Get-Date
$Duration2 = $End2 - $Start2
"Fin : $End2" | Out-File $LogFile -Append
"Durée : $($Duration2.ToString())" | Out-File $LogFile -Append
"" | Out-File $LogFile -Append

Write-Host "uk-2002 terminé en $($Duration2.TotalHours.ToString('F1')) heures" -ForegroundColor Yellow
Write-Host ""

# --- RESTAURATION SYSTÈME ---
Write-Host "=== Restauration du système ===" -ForegroundColor Cyan
powercfg /setactive SCHEME_BALANCED
powercfg /change standby-timeout-ac 30
powercfg /change monitor-timeout-ac 15
Write-Host "[OK] Paramètres d'alimentation restaurés"

# --- RÉSUMÉ FINAL ---
$TotalDuration = $End2 - $Start1  # Start1 = début Harvard500
Write-Host ""
Write-Host "=== Expériences terminées ===" -ForegroundColor Cyan
Write-Host "Durée totale : $($TotalDuration.TotalHours.ToString('F1')) heures"
Write-Host "Résultats dans : results.csv et classements.csv"
Write-Host "Log dans : $LogFile"

"=== RÉSUMÉ ===" | Out-File $LogFile -Append
"Durée totale : $($TotalDuration.ToString())" | Out-File $LogFile -Append
"Fin : $(Get-Date)" | Out-File $LogFile -Append
