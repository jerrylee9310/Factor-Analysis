import numpy as np
import pandas as pd
from typing import Optional, Literal, Dict, List, Union
from tqdm import tqdm
import matplotlib.pyplot as plt

class GenotypeValidator:
    """Genotype 데이터 검증을 담당하는 클래스"""
    
    @staticmethod
    def validate_maf(
        genotypes: np.ndarray,
        maf: np.ndarray,
        role: Optional[Literal["father", "mother", "offspring"]] = None,
        plot: bool = False
    ) -> Dict[str, Union[float, np.ndarray]]:
        """생성된 genotype의 MAF를 계산하여 초기 설정된 MAF와 비교

        Args:
            genotypes (np.ndarray): 검증할 genotype 데이터
            maf (np.ndarray): 원본 MAF 값
            role (Optional[Literal["father", "mother", "offspring"]]): 
                특정 역할의 genotype만 사용하여 검증. None일 경우 모든 구성원의 genotype 사용
            plot (bool): 결과를 시각화할지 여부. 기본값은 False

        Returns:
            Dict[str, Union[float, np.ndarray]]: 검증 결과를 담은 딕셔너리
                - correlation: 설정 MAF와 계산된 MAF 간의 상관계수
                - mean_diff: 설정 MAF와 계산된 MAF 간의 평균 절대 차이
                - original_maf: 원래 설정된 MAF 값들
                - computed_maf: genotype으로부터 계산된 MAF 값들
        """
        # 각 variant별 MAF 계산
        n_samples = len(genotypes)
        allele_count = np.sum(genotypes == 1, axis=0) + 2 * np.sum(genotypes == 2, axis=0)
        computed_maf = allele_count / (2 * n_samples)
        
        # 원본 MAF와 비교
        correlation = np.corrcoef(maf, computed_maf)[0, 1]
        mean_diff = np.mean(np.abs(maf - computed_maf))
        
        # 결과 출력
        role_text = f" ({role})" if role else ""
        print(f"MAF Validation Results{role_text}:")
        print(f"- Correlation: {correlation:.3f}")
        print(f"- Mean Absolute Difference: {mean_diff:.3f}")
        
        # plot 옵션이 True인 경우 시각화
        if plot:
            plt.figure(figsize=(10, 6))
            plt.scatter(maf, computed_maf, alpha=0.5)
            
            # 이상적인 라인 (x=y) 추가
            min_maf = min(maf.min(), computed_maf.min())
            max_maf = max(maf.max(), computed_maf.max())
            plt.plot([min_maf, max_maf], [min_maf, max_maf], 
                    'r--', label='Ideal (x=y)')
            
            plt.xlabel('Original MAF')
            plt.ylabel('Computed MAF')
            title = f'MAF Comparison{role_text}\nCorrelation: {correlation:.3f}, Mean Diff: {mean_diff:.3f}'
            plt.title(title)
            
            plt.legend()
            plt.grid(True)
            plt.show()
        
        return {
            "correlation": correlation,
            "mean_diff": mean_diff,
            "original_maf": maf,
            "computed_maf": computed_maf
        }
    
    @staticmethod
    def validate_correlation(
        geno_df: pd.DataFrame,
        plot: bool = False
    ) -> Dict[str, List[float]]:
        """가족 구성원 간의 genotype correlation을 계산하고 검증

        Args:
            geno_df (pd.DataFrame): 검증할 genotype 데이터
            plot (bool): 결과를 시각화할지 여부. 기본값은 False

        Returns:
            Dict[str, List[float]]: 검증 결과를 담은 딕셔너리
                - parent_correlation: 부모 간 correlation
                - father_offspring_correlations: 아버지-자녀 간 correlation 리스트
                - mother_offspring_correlations: 어머니-자녀 간 correlation 리스트
                - sibling_correlations: 자녀들 간의 correlation 리스트
        """
        results = {
            "parent_correlation": [],
            "father_offspring_correlations": [],
            "mother_offspring_correlations": [],
            "sibling_correlations": []
        }
        
        # 각 가족별로 correlation 계산
        family_ids = geno_df["family_id"].unique()
        for family_id in tqdm(family_ids, desc="Calculating family correlations"):
            family_data = geno_df[geno_df["family_id"] == family_id]
            
            # genotype 데이터 추출 (variant 열만 선택)
            variant_cols = [col for col in family_data.columns 
                           if col not in ["family_id", "role"]]
            
            # 부모 데이터 추출
            father_geno = family_data[family_data["role"] == "father"][variant_cols].values
            mother_geno = family_data[family_data["role"] == "mother"][variant_cols].values
            offspring_geno = family_data[family_data["role"] == "offspring"][variant_cols].values
            
            # 1. Parent correlation
            parent_corr = np.corrcoef(father_geno.flatten(), mother_geno.flatten())[0, 1]
            results["parent_correlation"].append(parent_corr)
            
            # 2. Parent-offspring correlations
            for off_idx in range(len(offspring_geno)):
                # Father-offspring
                father_off_corr = np.corrcoef(father_geno.flatten(), 
                                            offspring_geno[off_idx].flatten())[0, 1]
                results["father_offspring_correlations"].append(father_off_corr)
                
                # Mother-offspring
                mother_off_corr = np.corrcoef(mother_geno.flatten(), 
                                            offspring_geno[off_idx].flatten())[0, 1]
                results["mother_offspring_correlations"].append(mother_off_corr)
            
            # 3. Sibling correlations
            if len(offspring_geno) > 1:
                for i in range(len(offspring_geno)):
                    for j in range(i + 1, len(offspring_geno)):
                        sib_corr = np.corrcoef(offspring_geno[i].flatten(), 
                                             offspring_geno[j].flatten())[0, 1]
                        results["sibling_correlations"].append(sib_corr)
        
        # 결과 출력
        print("\nGenotype Correlation Validation Results:")
        print(f"Parent Correlation: {np.mean(results['parent_correlation']):.3f}")
        print(f"Father-Offspring Correlation: {np.mean(results['father_offspring_correlations']):.3f}")
        print(f"Mother-Offspring Correlation: {np.mean(results['mother_offspring_correlations']):.3f}")
        if results["sibling_correlations"]:
            print(f"Sibling Correlation: {np.mean(results['sibling_correlations']):.3f}")
        
        if plot:
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle('Genotype Correlations Distribution')
            
            # Parent correlation
            axes[0, 0].hist(results['parent_correlation'], bins=20)
            axes[0, 0].set_title('Parent Correlation')
            axes[0, 0].set_xlabel('Correlation')
            axes[0, 0].set_ylabel('Count')
            
            # Father-offspring correlation
            axes[0, 1].hist(results['father_offspring_correlations'], bins=20)
            axes[0, 1].set_title('Father-Offspring Correlation')
            axes[0, 1].set_xlabel('Correlation')
            
            # Mother-offspring correlation
            axes[1, 0].hist(results['mother_offspring_correlations'], bins=20)
            axes[1, 0].set_title('Mother-Offspring Correlation')
            axes[1, 0].set_xlabel('Correlation')
            axes[1, 0].set_ylabel('Count')
            
            # Sibling correlation
            if results["sibling_correlations"]:
                axes[1, 1].hist(results['sibling_correlations'], bins=20)
                axes[1, 1].set_title('Sibling Correlation')
                axes[1, 1].set_xlabel('Correlation')
            
            plt.tight_layout()
            plt.show()
        
        return results 