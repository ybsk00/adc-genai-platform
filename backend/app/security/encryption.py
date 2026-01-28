"""
AES-256-GCM Encryption Service
Premium 사용자 SMILES 암호화 저장
"""
from typing import Dict, Tuple
import os
import base64
import hashlib
import logging

logger = logging.getLogger(__name__)


class EncryptionService:
    """
    AES-256-GCM 암호화 서비스

    유료 사용자의 SMILES 코드를 암호화하여 저장
    - AES-256-GCM: 인증된 암호화
    - 사용자별 AAD (Additional Authenticated Data)
    - IV (Initialization Vector) 랜덤 생성
    """

    def __init__(self):
        self.key = self._load_encryption_key()

    def _load_encryption_key(self) -> bytes:
        """암호화 키 로드"""
        key_str = os.environ.get("ENCRYPTION_KEY")

        if not key_str:
            logger.warning("ENCRYPTION_KEY not set, using default (INSECURE!)")
            # 개발용 기본 키 (프로덕션에서는 절대 사용 금지!)
            key_str = base64.b64encode(b"dev_key_32_bytes_long_for_aes256").decode()

        try:
            return base64.b64decode(key_str)
        except Exception as e:
            logger.error(f"Failed to decode encryption key: {e}")
            raise ValueError("Invalid ENCRYPTION_KEY format")

    def encrypt_smiles(self, smiles: str, user_id: str) -> Dict[str, bytes]:
        """
        SMILES 암호화

        Args:
            smiles: 원본 SMILES 문자열
            user_id: 사용자 ID (AAD로 사용)

        Returns:
            {
                "encrypted_smiles": 암호화된 데이터,
                "encryption_iv": IV,
                "encryption_tag": 인증 태그,
                "structure_hash": 원본 해시 (검증용)
            }
        """
        try:
            from cryptography.hazmat.primitives.ciphers.aead import AESGCM

            # Generate random IV
            iv = os.urandom(12)

            # AAD: user_id 바인딩
            aad = user_id.encode()

            # Encrypt
            aesgcm = AESGCM(self.key)
            ciphertext = aesgcm.encrypt(iv, smiles.encode(), aad)

            # ciphertext에 tag가 포함됨 (마지막 16바이트)
            return {
                "encrypted_smiles": ciphertext[:-16],
                "encryption_iv": iv,
                "encryption_tag": ciphertext[-16:],
                "structure_hash": self._hash_smiles(smiles)
            }

        except ImportError:
            logger.warning("cryptography not installed, using fallback encoding")
            return self._fallback_encode(smiles, user_id)

        except Exception as e:
            logger.error(f"Encryption error: {e}")
            raise

    def decrypt_smiles(
        self,
        encrypted_data: bytes,
        iv: bytes,
        tag: bytes,
        user_id: str
    ) -> str:
        """
        SMILES 복호화

        Args:
            encrypted_data: 암호화된 데이터
            iv: Initialization Vector
            tag: 인증 태그
            user_id: 사용자 ID (AAD 검증)

        Returns:
            복호화된 SMILES 문자열
        """
        try:
            from cryptography.hazmat.primitives.ciphers.aead import AESGCM

            aad = user_id.encode()
            ciphertext_with_tag = encrypted_data + tag

            aesgcm = AESGCM(self.key)
            plaintext = aesgcm.decrypt(iv, ciphertext_with_tag, aad)

            return plaintext.decode()

        except ImportError:
            logger.warning("cryptography not installed, using fallback decoding")
            return self._fallback_decode(encrypted_data, user_id)

        except Exception as e:
            logger.error(f"Decryption error: {e}")
            raise

    def _hash_smiles(self, smiles: str) -> str:
        """복호화 후 검증용 해시"""
        return hashlib.sha256(smiles.encode()).hexdigest()

    def verify_smiles(self, smiles: str, expected_hash: str) -> bool:
        """복호화된 SMILES 검증"""
        return self._hash_smiles(smiles) == expected_hash

    def _fallback_encode(self, smiles: str, user_id: str) -> Dict[str, bytes]:
        """Fallback: Base64 인코딩 (개발용)"""
        encoded = base64.b64encode(smiles.encode())
        return {
            "encrypted_smiles": encoded,
            "encryption_iv": b"fallback",
            "encryption_tag": b"fallback",
            "structure_hash": self._hash_smiles(smiles)
        }

    def _fallback_decode(self, encoded_data: bytes, user_id: str) -> str:
        """Fallback: Base64 디코딩 (개발용)"""
        return base64.b64decode(encoded_data).decode()


class EncryptedStructureManager:
    """암호화된 구조 관리자"""

    def __init__(self):
        self.encryption_service = EncryptionService()
        from app.core.supabase import get_supabase_client
        self.supabase = get_supabase_client()

    async def store_encrypted_structure(
        self,
        session_id: str,
        user_id: str,
        smiles: str,
        structure_type: str = "candidate",
        rank: int = 1
    ) -> str:
        """
        암호화된 구조 저장

        Returns:
            저장된 레코드 ID
        """
        encrypted = self.encryption_service.encrypt_smiles(smiles, user_id)

        result = await self.supabase.table("encrypted_structures").insert({
            "session_id": session_id,
            "user_id": user_id,
            "encrypted_smiles": encrypted["encrypted_smiles"],
            "encryption_iv": encrypted["encryption_iv"],
            "encryption_tag": encrypted["encryption_tag"],
            "structure_type": structure_type,
            "rank": rank,
            "structure_hash": encrypted["structure_hash"]
        }).execute()

        return result.data[0]["id"] if result.data else None

    async def retrieve_decrypted_structure(
        self,
        structure_id: str,
        user_id: str
    ) -> str:
        """
        암호화된 구조 조회 및 복호화

        Args:
            structure_id: 구조 ID
            user_id: 사용자 ID (권한 및 AAD 검증)

        Returns:
            복호화된 SMILES
        """
        result = await self.supabase.table("encrypted_structures").select("*").eq(
            "id", structure_id
        ).eq(
            "user_id", user_id  # 소유권 확인
        ).single().execute()

        if not result.data:
            raise ValueError("Structure not found or access denied")

        data = result.data

        # 복호화
        smiles = self.encryption_service.decrypt_smiles(
            encrypted_data=data["encrypted_smiles"],
            iv=data["encryption_iv"],
            tag=data["encryption_tag"],
            user_id=user_id
        )

        # 해시 검증
        if not self.encryption_service.verify_smiles(smiles, data["structure_hash"]):
            raise ValueError("Structure integrity check failed")

        # 접근 시간 업데이트
        await self.supabase.table("encrypted_structures").update({
            "accessed_at": "now()"
        }).eq("id", structure_id).execute()

        return smiles
