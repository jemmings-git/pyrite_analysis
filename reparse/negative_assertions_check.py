import spacy
import pandas

NLP = spacy.load('en')


def check_assertions(phrases):
    for phrase in phrases:
        if not 'not' in phrase:
            continue
        negative_chain(phrase)


def negative_chain(phrase):
    doc = NLP(phrase)
    for token in doc:
        arc = [child.text.lower() for child in token.children]
        if 'not' in arc and 'pyrite' in arc:
            print(doc)
            print(token.text, token.dep_, token.head.text, token.head.pos_,
                 [child for child in token.children])


if __name__ == '__main__':
    results = pandas.read_csv('../data/results.csv')
    check_assertions(results['phrase'])
