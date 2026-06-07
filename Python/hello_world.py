# Hello world script: learn using LangGraph to create agentic workflows
# Reva S
# 06-Jun-2026
# -----------

# load requirements
from langgraph.graph import StateGraph, MessagesState, START, END

# define function
def mock_llm(state: MessagesState):
    return {"messages": [{"role": "ai", "content": "hello world"}]}

# create data
graph = StateGraph(MessagesState)
graph.add_node(mock_llm)
graph.add_edge(START, "mock_llm")
graph.add_edge("mock_llm", END)
graph = graph.compile()

# run function
response = graph.invoke({"messages": [{"role": "user", "content": "hi!"}]})
# output
print(response)