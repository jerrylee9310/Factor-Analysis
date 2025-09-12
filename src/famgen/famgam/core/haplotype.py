from dataclasses import dataclass, field
import numpy as np
from typing import Optional, Tuple

@dataclass
class HaplotypeGenerator:
    """Haplotype 생성 및 관리를 담당하는 클래스
    
    Attributes:
        num_variant (int): 생성할 변이의 수
        maf_lim (list): MAF(Minor Allele Frequency) 범위 [min, max]
    """
    num_variant: int
    maf_lim: list = field(default_factory=lambda: [0.05, 0.95])
    
    # 내부 변수
    maf: np.ndarray = field(init=False)
    father_haplotype: Optional[np.ndarray] = field(init=False, default=None)
    mother_haplotype: Optional[np.ndarray] = field(init=False, default=None)
    offspring_haplotype: Optional[np.ndarray] = field(init=False, default=None)

    def __post_init__(self):
        """초기화 후 MAF 생성"""
        self.maf = np.random.uniform(
            min(self.maf_lim), 
            max(self.maf_lim), 
            self.num_variant
        )
    
    def __generate_haplotype(self, num_sample: int) -> np.ndarray:
        """단일 샘플의 haplotype 생성
        
        Args:
            num_sample (int): 생성할 샘플 수
            
        Returns:
            np.ndarray: (num_sample, num_variant, 2) 형태의 haplotype 행렬
        """
        n = num_sample
        m = len(self.maf)
        
        hap1 = np.random.binomial(n=1, p=self.maf, size=(n, m)).astype(int)
        hap2 = np.random.binomial(n=1, p=self.maf, size=(n, m)).astype(int)
        return np.dstack([hap1, hap2])

    def __random_inherit(self, haps_parent: np.ndarray) -> np.ndarray:
        """부모의 haplotype에서 무작위로 하나를 상속
        
        Args:
            haps_parent (np.ndarray): 부모의 haplotype
            
        Returns:
            np.ndarray: 상속된 haplotype
        """
        n, m, _ = haps_parent.shape
        # 각 위치에서 첫 번째(0) 또는 두 번째(1) haplotype을 선택
        mask = np.random.randint(2, size=(n, m))
        # 선택된 haplotype 값을 가져옴
        return np.take_along_axis(haps_parent, mask[..., np.newaxis], axis=2).squeeze(axis=2)
    
    def generate_parent_haplotype(self, num_sample: int, return_haplotypes: bool = False) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """부모의 haplotype 생성 및 저장
        
        Args:
            num_sample (int): 생성할 샘플 수
            return_haplotypes (bool): haplotype을 반환할지 여부. 기본값은 False
            
        Returns:
            Optional[Tuple[np.ndarray, np.ndarray]]: 
                return_haplotypes가 True일 경우 (father_haplotype, mother_haplotype) 반환
                False일 경우 None 반환
        """
        self.father_haplotype = self.__generate_haplotype(num_sample)
        self.mother_haplotype = self.__generate_haplotype(num_sample)
        
        if return_haplotypes:
            return self.father_haplotype, self.mother_haplotype
        return None
    
    def get_father_haplotype(self) -> np.ndarray:
        """아버지의 haplotype 반환
        
        Returns:
            np.ndarray: 아버지의 haplotype
        """
        if self.father_haplotype is None:
            raise ValueError("부모의 haplotype이 아직 생성되지 않았습니다. generate_parent_haplotype을 먼저 실행해주세요.")
        return self.father_haplotype
    
    def get_mother_haplotype(self) -> np.ndarray:
        """어머니의 haplotype 반환
        
        Returns:
            np.ndarray: 어머니의 haplotype
        """
        if self.mother_haplotype is None:
            raise ValueError("부모의 haplotype이 아직 생성되지 않았습니다. generate_parent_haplotype을 먼저 실행해주세요.")
        return self.mother_haplotype
    
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
                
        Raises:
            ValueError: 부모의 haplotype이 주어지지 않고 내부에도 저장되어 있지 않은 경우
        """
        # 부모 haplotype 결정
        if haps_mother is None:
            if self.mother_haplotype is None:
                raise ValueError("어머니의 haplotype이 주어지지 않았고 내부에도 저장되어 있지 않습니다.")
            haps_mother = self.mother_haplotype
            
        if haps_father is None:
            if self.father_haplotype is None:
                raise ValueError("아버지의 haplotype이 주어지지 않았고 내부에도 저장되어 있지 않습니다.")
            haps_father = self.father_haplotype
            
        # 입력값 검증
        if haps_mother.shape != haps_father.shape:
            raise ValueError("부모의 haplotype shape이 서로 다릅니다.")
            
        n_parents, n_variants, _ = haps_mother.shape
        
        # 각 부모 쌍에 대해 자손 생성
        offspring_list = []
        for i in range(n_parents):
            parent_mother = haps_mother[i:i+1]  # (1, n_variants, 2) 형태 유지
            parent_father = haps_father[i:i+1]
            
            # num_offspring만큼 자손 생성
            for _ in range(num_offspring):
                # 어머니로부터 상속
                hap_from_mom = self.__random_inherit(parent_mother)
                # 아버지로부터 상속
                hap_from_dad = self.__random_inherit(parent_father)
                
                # (1, n_variants, 2) 형태로 만들어서 저장
                offspring = np.dstack([hap_from_mom, hap_from_dad])
                offspring_list.append(offspring)
        
        # 모든 자손을 하나의 배열로 합침
        self.offspring_haplotype = np.vstack(offspring_list)
        
        if return_haplotypes:
            return self.offspring_haplotype
        return None

    def get_offspring_haplotype(self) -> np.ndarray:
        """생성된 자손의 haplotype 반환
            
        Returns:
            np.ndarray: 자손의 haplotype
            
        Raises:
            ValueError: 자손의 haplotype이 아직 생성되지 않은 경우
        """
        if self.offspring_haplotype is None:
            raise ValueError("자손의 haplotype이 아직 생성되지 않았습니다.")
        return self.offspring_haplotype 