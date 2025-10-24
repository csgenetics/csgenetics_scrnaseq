from .unified_qc import MultiqcModule
import os
import yaml
from multiqc import config

def register_plugin_search_patterns():
    """Hook called by MultiQC before config is loaded.

    This registers our search patterns so MultiQC knows how to find our files.
    """
    # Load our config file from the unified_qc package directory
    config_file = os.path.join(os.path.dirname(__file__), 'multiqc_config.yaml')
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            plugin_config = yaml.safe_load(f)
            if plugin_config and 'sp' in plugin_config:
                # Merge our search patterns into MultiQC's config
                if not hasattr(config, 'sp') or config.sp is None:
                    config.sp = {}
                config.sp.update(plugin_config['sp'])
