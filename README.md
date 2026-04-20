# ADAR-sensor: AI-powered RNA Sensor Design Platform

ViennaRNA 물리 엔진과 AI 추론을 결합한 ADAR 편집 센서 설계 및 최적화 플랫폼입니다.

## 주요 구성 요소

### 1. ViennaRNA 시각화 스킬 (`.gemini/skills/viennarna-visualization`)
- **핵심 기능**: 거대한 mRNA 내에서 센서가 결합하는 dsRNA 하이브리드 구조 시각화
- **작동 원리**:
  1. `duplexfold`로 최적 결합 위치 탐색 (MFE 기반)
  2. 점-괄호(Dot-bracket) 표기법 해석
  3. 커스텀 matplotlib 렌더링

### 2. C3 Sensor 프로젝트 (`Projects/C3_Sensor/`)
- ALU 기반 센서 설계 및 검증
- dsRNA 이량체 시각화
-候选자 점수화 및 순위 결정

## 디렉토리 구조

```
ADAR-sensor/
├── SKILL.md                 # ViennaRNA 시각화 스킬 설명
├── manifest.json            # 스킬 메타데이터
├── scripts/
│   ├── visualize_binding.py # 핵심 시각화 엔진
│   ├── viennarna_enhanced.py
│   └── plot_native_dimer.py
└── c3_sensor/
    ├── plot_dsRNA.py
    ├── plot_native_dimer.py
    ├── plot_native_dimer_styled.py
    ├── visualize_c3_elite.py
    ├── c3_300bp_candidates_scored.csv
    ├── top6_c3_sites.txt
    └── human_C3_mRNA.fasta
```

## 사용법

### 센서 결합 시각화
```bash
python scripts/visualize_binding.py --mRNA <target_mRNA.fasta> --sensor <sensor_sequence>
```

### C3 센서 분석
```bash
python c3_sensor/visualize_c3_elite.py
```

## 과학적 배경

- **ADAR (Adenosine Deaminase Acting on RNA)**: 전사후 RNA 편집 효소
- **C3 ALU 미메틱 센서**: ALU 서열 기반 설계된 센서
- **ViennaRNA**: RNA 이차구조 예측 및 분석을 위한 물리화학적 모델

## 라이선스

연구용으로만 사용하세요.