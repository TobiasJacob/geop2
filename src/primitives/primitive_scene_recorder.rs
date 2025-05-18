use crate::algebra_error::AlgebraResult;

use super::primitive_scene::PrimitiveScene;

pub struct PrimitiveSceneRecorder {
    pub scenes: Vec<PrimitiveScene>,
}

impl PrimitiveSceneRecorder {
    pub fn new() -> Self {
        Self { scenes: vec![] }
    }

    pub fn add_scene(&mut self, scene: PrimitiveScene) {
        self.scenes.push(scene);
    }

    pub fn save_to_folder(&self, folder_path: &str) -> AlgebraResult<()> {
        for (i, scene) in self.scenes.iter().enumerate() {
            let file_path = format!("{}/scene_{}.html", folder_path, i);
            scene.save_to_file(&file_path)?;
        }
        println!("Saved {} scenes to {}", self.scenes.len(), folder_path);
        Ok(())
    }
}
