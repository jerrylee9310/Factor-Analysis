# FAMGAM

FAMGAM은 부모의 유전형을 기반으로 자손의 유전형을 생성하는 Python 패키지입니다.

## 설치 방법

```bash
pip install famgam
```

## 사용 예시

```python
from famgam.core import FamGenoSimul, GenotypeValidator

# 시뮬레이터 초기화
fam_geno = FamGenoSimul(num_variant=1000)

# 부모 haplotype 생성
fam_geno.generate_parent_haplotype(num_sample=100)

# 자손 haplotype 생성
fam_geno.make_offspring_haplotype(num_offspring=2)

# Genotype 데이터 가져오기
geno_df = fam_geno.get_genotype(as_dataframe=True, std=True)

# MAF 검증
validator = GenotypeValidator()
maf_results = validator.validate_maf(
    genotypes=geno_df.drop(["family_id", "role"], axis=1).values,
    maf=fam_geno.haplotype_generator.maf,
    plot=True
)

# 상관관계 검증
corr_results = validator.validate_correlation(geno_df, plot=True)
```

## 주요 기능

### FamGenoSimul

- `__init__(num_variant: int, maf_lim: Optional[List[float]] = None)`: 시뮬레이터 초기화
- `generate_parent_haplotype(num_sample: int, return_haplotypes: bool = False)`: 부모의 haplotype 생성
- `make_offspring_haplotype(num_offspring: int = 1, haps_mother: Optional[np.ndarray] = None, haps_father: Optional[np.ndarray] = None, return_haplotypes: bool = False)`: 자손의 haplotype 생성
- `get_genotype(role: Optional[Literal["father", "mother", "offspring"]] = None, as_dataframe: bool = True, std: bool = False)`: Genotype 데이터 가져오기

### GenotypeValidator

- `validate_maf(genotypes: np.ndarray, maf: np.ndarray, role: Optional[Literal["father", "mother", "offspring"]] = None, plot: bool = False)`: MAF 검증
- `validate_correlation(geno_df: pd.DataFrame, plot: bool = False)`: 상관관계 검증

## 요구사항

- Python >= 3.7
- NumPy >= 1.20.0
- Pandas >= 1.3.0
- Matplotlib >= 3.4.0
- tqdm >= 4.62.0 