"""
Prepare gene embeddings for the training pipeline.
"""


class BuildTrain:
    def __init__(self, mode="load", model=None):

        if mode not in ["load", "generate"]:
            raise ValueError("Mode must be 'load' or 'generate'")

        self.mode = mode

        self.embeddings = None
        self.gene_descriptions = None
        self.model = model

    def load_precomputed(self, embeddings):
        """Use embeddings already loaded by LoadTrain."""
        self.embeddings = embeddings
        return self

    def load_gene_descriptions(self, gene_descriptions):
        """Use gene descriptions already loaded by LoadTrain."""
        self.gene_descriptions = gene_descriptions
        return self

    def load_model(self):
        """Load or initialize the embedding model."""
        # To be implemented once the embedding model is selected.
        raise NotImplementedError

    def generate_embeddings(self):
        """Generate embeddings from gene descriptions."""
        # To be implemented once the embedding model is selected.
        raise NotImplementedError

    def validate_embeddings(self):
        """Validate the structure and dimensions of the embeddings."""
        # To be implemented once the embedding model is selected.
        raise NotImplementedError

    def build(self):

        if self.mode == "load":
            self.load_precomputed()

        elif self.mode == "generate":
            self.load_gene_descriptions()

            if self.model is None:
                self.load_model()

            self.generate_embeddings()

        self.validate_embeddings()

        return self
