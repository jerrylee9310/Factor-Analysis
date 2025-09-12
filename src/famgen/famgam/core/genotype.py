import numpy as np
import pandas as pd
from typing import Union, Optional, Literal, List
from .haplotype import HaplotypeGenerator

class GenotypeConverter:
    """Haplotype을 Genotype으로 변환하는 클래스"""
    
    @staticmethod
    def haplotype_to_genotype(haplotype: np.ndarray) -> np.ndarray:
        """Haplotype을 Genotype으로 변환
        
        Args:
            haplotype (np.ndarray): (n_samples, n_variants, 2) 형태의 haplotype
            
        Returns:
            np.ndarray: (n_samples, n_variants) 형태의 genotype
        """
        return np.sum(haplotype, axis=2)
    
    @staticmethod
    def standardize_genotype(genotype: np.ndarray, maf: np.ndarray) -> np.ndarray:
        """Genotype을 MAF 기준으로 standardize

        Args:
            genotype (np.ndarray): standardize할 genotype 데이터
            maf (np.ndarray): 각 variant의 MAF 값

        Returns:
            np.ndarray: standardized genotype 데이터
        """
        # Expected mean: 2 * MAF (각 variant position별)
        expected_mean = 2 * maf
        
        # Expected variance: 2 * MAF * (1 - MAF)
        expected_std = np.sqrt(2 * maf * (1 - maf))
        
        # Standardization
        standardized = (genotype - expected_mean) / expected_std
        
        return standardized

class FamGenoSimul:
    """가족 구성원의 유전형을 생성하는 클래스
    
    Attributes:
        num_variant (int): 생성할 변이의 수
        maf_lim (list): MAF(Minor Allele Frequency) 범위 [min, max]
    """
    def __init__(self, num_variant: int, maf_lim: Optional[List[float]] = None):
        if maf_lim is None:
            maf_lim = [0.05, 0.95]
            
        self.haplotype_generator = HaplotypeGenerator(num_variant, maf_lim)
        self.genotype_converter = GenotypeConverter()
    
    def generate_parent_haplotype(self, num_sample: int, return_haplotypes: bool = False) -> Optional[tuple[np.ndarray, np.ndarray]]:
        """부모의 haplotype 생성 및 저장
        
        Args:
            num_sample (int): 생성할 샘플 수
            return_haplotypes (bool): haplotype을 반환할지 여부. 기본값은 False
            
        Returns:
            Optional[tuple[np.ndarray, np.ndarray]]: 
                return_haplotypes가 True일 경우 (father_haplotype, mother_haplotype) 반환
                False일 경우 None 반환
        """
        return self.haplotype_generator.generate_parent_haplotype(num_sample, return_haplotypes)
    
    def make_offspring_haplotype(
        self, 
        num_offspring: int = 1,
        haps_mother: Optional[np.ndarray] = None, 
        haps_father: Optional[np.ndarray] = None,
        return_haplotypes: bool = False
    ) -> Optional[np.ndarray]:
        """자손의 haplotype 생성
        
        Args:
            num_offspring (int): 각 부모 쌍당 생성할 자손의 수. 기본값은 1
            haps_mother (Optional[np.ndarray]): 어머니의 haplotype. 
                None일 경우 내부 저장된 값 사용
            haps_father (Optional[np.ndarray]): 아버지의 haplotype. 
                None일 경우 내부 저장된 값 사용
            return_haplotypes (bool): haplotype을 반환할지 여부. 기본값은 False
            
        Returns:
            Optional[np.ndarray]: 
                return_haplotypes가 True일 경우 offspring_haplotype 반환
                False일 경우 None 반환
        """
        return self.haplotype_generator.make_offspring_haplotype(
            num_offspring, haps_mother, haps_father, return_haplotypes
        )
    
    def get_maf(self) -> np.ndarray:
        """MAF(Minor Allele Frequency) 값을 반환
        
        Returns:
            np.ndarray: 각 variant의 MAF 값
        """
        return self.haplotype_generator.maf
    
    def get_genotype(self, 
                    role: Optional[Literal["father", "mother", "offspring"]] = None,
                    as_dataframe: bool = True,
                    std: bool = False
                    ) -> Union[np.ndarray, pd.DataFrame]:
        """Haplotype을 Genotype으로 변환하여 반환
        
        Args:
            role (Optional[Literal["father", "mother", "offspring"]]): 
                특정 역할의 genotype만 반환. None일 경우 모든 구성원의 genotype 반환
            as_dataframe (bool): DataFrame 형식으로 반환할지 여부. 기본값은 True
            std (bool): Standardized genotype을 반환할지 여부. 기본값은 False
            
        Returns:
            Union[np.ndarray, pd.DataFrame]: 
                as_dataframe이 True이면 DataFrame, False이면 numpy array 반환
                
        Raises:
            ValueError: 요청된 역할의 haplotype이 아직 생성되지 않은 경우
        """
        genotypes = []
        roles = []
        family_ids = []
        
        # 아버지 genotype
        if (role is None or role == "father"):
            father_hap = self.haplotype_generator.get_father_haplotype()
            father_geno = self.genotype_converter.haplotype_to_genotype(father_hap)
            if std:
                father_geno = self.genotype_converter.standardize_genotype(
                    father_geno, self.haplotype_generator.maf
                )
            n_fathers = len(father_geno)
            genotypes.append(father_geno)
            roles.extend(["father"] * n_fathers)
            family_ids.extend(range(n_fathers))
            
        # 어머니 genotype
        if (role is None or role == "mother"):
            mother_hap = self.haplotype_generator.get_mother_haplotype()
            mother_geno = self.genotype_converter.haplotype_to_genotype(mother_hap)
            if std:
                mother_geno = self.genotype_converter.standardize_genotype(
                    mother_geno, self.haplotype_generator.maf
                )
            n_mothers = len(mother_geno)
            genotypes.append(mother_geno)
            roles.extend(["mother"] * n_mothers)
            family_ids.extend(range(n_mothers))
            
        # 자손 genotype
        if (role is None or role == "offspring"):
            offspring_hap = self.haplotype_generator.get_offspring_haplotype()
            offspring_geno = self.genotype_converter.haplotype_to_genotype(offspring_hap)
            if std:
                offspring_geno = self.genotype_converter.standardize_genotype(
                    offspring_geno, self.haplotype_generator.maf
                )
            n_offspring = len(offspring_geno)
            n_families = len(father_hap) if father_hap is not None else len(mother_hap)
            offspring_per_family = n_offspring // n_families
            
            genotypes.append(offspring_geno)
            roles.extend(["offspring"] * n_offspring)
            family_ids.extend([i // offspring_per_family for i in range(n_offspring)])
            
        if not genotypes:
            raise ValueError("요청된 역할의 haplotype이 아직 생성되지 않았습니다.")
            
        # 모든 genotype을 하나의 배열로 합침
        all_genotypes = np.vstack(genotypes)
        
        if not as_dataframe:
            return all_genotypes
            
        # DataFrame 생성
        variant_cols = [f"variant_{i+1}" for i in range(all_genotypes.shape[1])]
        df = pd.DataFrame(all_genotypes, columns=variant_cols)
        
        # annotation 추가
        df.insert(0, "family_id", family_ids)
        df.insert(1, "role", roles)
        
        return df 