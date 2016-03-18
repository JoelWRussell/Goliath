/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;


import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.*;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

public class XMLHelper {
private Document doc;
private XPath path;

public XMLHelper()
{


}

public void setDocument(String sz) throws Exception{

    byte[] data = null;
    ByteArrayOutputStream ba = new ByteArrayOutputStream();
    try (PrintWriter pw = new PrintWriter(ba)) {
        pw.write(sz);
        pw.flush();
        data = ba.toByteArray();
    }
    if (data != null){
        ByteArrayInputStream is = new ByteArrayInputStream(data);
        setDocument(is);

    }
}
public void setDocument(File f) throws Exception{
    //parse the XML file so that it can be consumed
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        //Document is the root element
        doc = dBuilder.parse(f);
        doc.getDocumentElement().normalize();

        XPathFactory xpfactory = XPathFactory.newInstance();
        path = xpfactory.newXPath();
}
public void setDocument(InputStream is) throws Exception{
    //parse the XML file so that it can be consumed
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        //Document is the root element
        doc = dBuilder.parse(is);
        doc.getDocumentElement().normalize();

        XPathFactory xpfactory = XPathFactory.newInstance();
        path = xpfactory.newXPath();
}

public String EvaluateXPath(String xpath) throws XPathExpressionException{

    return path.evaluate(xpath, doc);

}
public String EvaluateXPathAttr(String xpath, String attr, int ix) throws XPathExpressionException{

    NodeList nodes =  (NodeList)path.evaluate(xpath, doc, XPathConstants.NODESET);
    Node child = nodes.item(ix);
    if (child instanceof Element){
        Element el = (Element)child;
        return el.getAttribute(attr);
    }
    return null;

}


public int EvaluateXPathInt(String xpath) throws XPathExpressionException{

    return ((Number)path.evaluate(xpath, doc, XPathConstants.NUMBER)).intValue();

}
public double EvaluateXPathDouble(String xpath) throws XPathExpressionException{

    return ((Number)path.evaluate(xpath, doc, XPathConstants.NUMBER)).doubleValue();

}
long EvaluateXPathLong(String xpath) throws XPathExpressionException {
    return ((Number)path.evaluate(xpath, doc, XPathConstants.NUMBER)).longValue();
}

}



