from sklearn.linear_model import Ridge


class ModelFactory:
    @staticmethod
    def create_model(model_name, **kwargs):

        models = {"ridge": Ridge}

        if model_name not in models:
            raise ValueError(f"Unknown model: {model_name}")

        return models[model_name](**kwargs)
