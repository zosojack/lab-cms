#!/usr/bin/env python3
"""
Script per convertire tutte le immagini PNG in PDF
nella cartella images e sue sottocartelle
"""

import os
from PIL import Image
import glob

def convert_png_to_pdf(input_path, output_path=None):
    """
    Converte un file PNG in PDF
    """
    try:
        # Apri l'immagine PNG
        with Image.open(input_path) as img:
            # Converti in RGB se necessario (i PNG possono avere trasparenza)
            if img.mode in ('RGBA', 'LA', 'P'):
                # Crea uno sfondo bianco
                rgb_img = Image.new('RGB', img.size, (255, 255, 255))
                if img.mode == 'P':
                    img = img.convert('RGBA')
                rgb_img.paste(img, mask=img.split()[-1] if img.mode in ('RGBA', 'LA') else None)
                img = rgb_img
            elif img.mode != 'RGB':
                img = img.convert('RGB')
            
            # Se non specificato, crea il nome del file PDF
            if output_path is None:
                output_path = input_path.rsplit('.', 1)[0] + '.pdf'
            
            # Salva come PDF
            img.save(output_path, 'PDF', resolution=100.0)
            print(f"✅ Convertito: {input_path} → {output_path}")
            return True
            
    except Exception as e:
        print(f"❌ Errore convertendo {input_path}: {e}")
        return False

def convert_all_images_in_directory(base_path):
    """
    Trova e converte tutte le immagini PNG nella directory e sottodirectory
    """
    print(f"🔍 Cercando immagini PNG in: {base_path}")
    
    # Pattern per trovare tutti i file PNG ricorsivamente
    pattern = os.path.join(base_path, "**", "*.png")
    png_files = glob.glob(pattern, recursive=True)
    
    # Aggiungi anche i PNG nella directory principale
    pattern_root = os.path.join(base_path, "*.png")
    png_files.extend(glob.glob(pattern_root))
    
    if not png_files:
        print("❌ Nessun file PNG trovato!")
        return
    
    print(f"📁 Trovati {len(png_files)} file PNG da convertire:")
    
    success_count = 0
    for png_file in png_files:
        print(f"🔄 Convertendo: {os.path.basename(png_file)}")
        if convert_png_to_pdf(png_file):
            success_count += 1
    
    print(f"\n🎉 Conversione completata!")
    print(f"✅ Convertiti con successo: {success_count}/{len(png_files)} file")
    
    if success_count < len(png_files):
        print(f"⚠️  Falliti: {len(png_files) - success_count} file")

if __name__ == "__main__":
    # Percorso della cartella images
    images_dir = "/Users/zosojack/lab-cms/images"
    
    print("🖼️  CONVERTITORE PNG → PDF")
    print("=" * 40)
    
    # Verifica che la directory esista
    if not os.path.exists(images_dir):
        print(f"❌ Directory non trovata: {images_dir}")
        exit(1)
    
    # Converti tutte le immagini
    convert_all_images_in_directory(images_dir)
