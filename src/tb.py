import tensorflow as tf
from tensorflow.core.framework import graph_pb2

# Load the graph
with tf.io.gfile.GFile('exported_graph.pb', 'rb') as f:
    graph_def = graph_pb2.GraphDef()
    graph_def.ParseFromString(f.read())

# Import the graph to a new Graph
with tf.Graph().as_default() as graph:
    tf.import_graph_def(graph_def, name="")

# Use TensorFlow 2.x summary writer
with tf.summary.create_file_writer('logs').as_default():
    tf.summary.graph(graph)
    tf.summary.flush()

# Now, run TensorBoard pointing it to the 'logs' directory.
